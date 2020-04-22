/* Code converted from shallow_base.f90 using F2C-ACC program. 
 * Manually replaced: 
 * - WRITE statements with printf
 * - MOD operator with % 
 * - system_clock with wtime
 * Fixed several of the array references which had x dimension as 1, 
 * instead of M_LEN. 
 * Fixed values set using d and e notation. 
 * (7 June 2011)
 ***************
 * 'Pure' C version developed by G.D Riley (UoM) (25 Jan 2012)
 * removed all ftocmacros
 * used sin and cos not sinf and cosf (since all data are doubles)
 * needed to declare arrays +1 to cope with Fortran indexing
 * Compile, e.g.:
 * gcc -O2 -fopenmp -o sb shallow_base_openmp_v3.c wtime.c -lm
 ** NOTE:  May need to set 'ulimit -s unlimited' to run large
 * problems (e.g. 512x512).
 * Results are consistent with Fortran version of the code.
 * GDR: July 2013
 * Applied static, chunk scheduling for load balance,as described
 * in Michail Pappas' MSc thesis (2012).
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
    #include <omp.h>
#endif

#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

#define TRUE 1
#define FALSE 0
#define M 128
#define N 128
#define M_LEN M + 1
#define N_LEN N + 1
#define ITMAX 2000
#define L_OUT TRUE

extern double wtime(); 

//===================================================

double compute_checksum(double field[M_LEN][N_LEN], 
			int lenx, int leny)
{
  int i, j;
  double sum = 0.0;

  for(i=0;i<lenx;i++){
    for(j=0;j<leny;j++){
      sum += field[i][j];
    }
  }

  return sum;
}

//===================================================

//! Benchmark weather prediction program for comparing the
//! preformance of current supercomputers. The model is
//! based on the paper - The Dynamics of Finite-Difference
//! Models of the Shallow-Water Equations, by Robert Sadourny
//! J. Atm. Sciences, Vol 32, No 4, April 1975.
//!     
//! Code by Paul N. Swarztrauber, National Center for
//! Atmospheric Research, Boulder, Co,  October 1984.
//! Modified by Juliana Rew, NCAR, January 2006
//!
//! In this version, shallow4.f, initial and calculated values
//! of U, V, and P are written to a netCDF file
//! for later use in visualizing the results. The netCDF data
//! management library is freely available from
//! http://www.unidata.ucar.edu/software/netcdf
//! This code is still serial but has been brought up to modern
//! Fortran constructs and uses portable intrinsic Fortran 90 timing routines. 
//! This can be compiled on the IBM SP using:
//! xlf90 -qmaxmem=-1 -g -o shallow4 -qfixed=132 -qsclk=micro \
//! -I/usr/local/include shallow4.f -L/usr/local/lib32/r4i4 -l netcdf
//! where the -L and -I point to local installation of netCDF
//!     
//! Changes from shallow4.f (Annette Osprey, January 2010):
//! - Converted to free-form fortran 90.  
//! - Some tidying up of old commented-out code.   
//! - Explicit type declarations.
//! - Variables n, m, ITMAX and mprint read in from namelist. 
//! - Dynamic array allocation.
//! - Only write to netcdf at mprint timesteps.
//! - Don't write wrap-around points to NetCDF file.
//! - Use 8-byte reals.
//!
//! Further changes (Annette Osprey & Graham Riley, February 2011): 
//! - Remove unnecessary halo updates.
//! - Split loops to improve TLB access.
//! - Modify timers to accumulate loop times over whole run. 
//! - Remove old-style indentation. 
//!
//! Minimal serial version (26 May 2011)

int main(int argc, char **argv) {
  
  // solution arrays
  double u[M_LEN][N_LEN],v[M_LEN][N_LEN],p[M_LEN][N_LEN];
  double unew[M_LEN][N_LEN],vnew[M_LEN][N_LEN],pnew[M_LEN][N_LEN];
  double uold[M_LEN][N_LEN],vold[M_LEN][N_LEN],pold[M_LEN][N_LEN];
  double cu[M_LEN][N_LEN],cv[M_LEN][N_LEN],z[M_LEN][N_LEN],h[M_LEN][N_LEN],psi[M_LEN][N_LEN];

  double dt,tdt,dx,dy,a,alpha,el,pi;
  double tpi,di,dj,pcf;
  double tdts8,tdtsdx,tdtsdy,fsdx,fsdy;

  int mnmin,ncycle;
  int i,j;
  int nthreads, chunk_size;

  // timer variables 
  double t100,t200,t300;
  double tstart,ctime,tcyc,time;
  double t100i,t200i,t300i;
  double c1,c2;

  // ** Initialisations ** 

  nthreads = 1;
#ifdef _OPENMP
  nthreads = omp_get_max_threads();
#endif
  chunk_size = (int) ceil( (float)M / (float)nthreads);

  // Note below that two delta t (tdt) is set to dt on the first
  // cycle after which it is reset to dt+dt.
  dt = 90.;
  tdt = dt;
 
  dx = 100000.;
  dy = 100000.;
  fsdx = 4. / dx;
  fsdy = 4. / dy;

  a = 1000000.;
  alpha = .001;

  el = N * dx;
  pi = 4. * atanf(1.);
  tpi = pi + pi;
  di = tpi / M;
  dj = tpi / N;
  pcf = pi * pi * a * a / (el * el);

  // Initial values of the stream function and p
#pragma omp parallel for default (shared) private(i,j)
  for (i=0;i<M_LEN;i++) {
    for (j=0;j<N_LEN;j++) {
      psi[i][j] = a * sin((i + .5) * di) * sin((j + .5) * dj);
      p[i][j] = pcf * (cos(2. * (i) * di) + cos(2. * (j) * dj)) + 50000.;
    }
  }
    
  // Initialize velocities
  #pragma omp parallel for default (shared) private(i,j)
  for (i=0;i<M;i++) {
    for (j=0;j<N;j++) {
      u[i + 1][j] = -(psi[i + 1][j + 1] - psi[i + 1][j]) / dy;
      v[i][j + 1] = (psi[i + 1][j + 1] - psi[i][j + 1]) / dx;
    }
  }
     
  // Periodic continuation
  for (j=0;j<N;j++) {
    u[0][j] = u[M][j];
    v[M][j + 1] = v[0][j + 1];
  }
  for (i=0;i<M;i++) {
    u[i + 1][N] = u[i + 1][0];
    v[i][0] = v[i][N];
  }
  u[0][N] = u[M][0];
  v[M][0] = v[0][N];
  #pragma omp parallel default (shared) private(i,j)
  for (i=0;i<M_LEN;i++) {
    for (j=0;j<N_LEN;j++) {
      uold[i][j] = u[i][j];
      vold[i][j] = v[i][j];
      pold[i][j] = p[i][j];
    }
  }
     
  // Print initial values
  if ( L_OUT ) {
    printf(" number of points in the x direction %d\n", N); 
    printf(" number of points in the y direction %d\n", M); 
    printf(" grid spacing in the x direction     %f\n", dx); 
    printf(" grid spacing in the y direction     %f\n", dy); 
    printf(" time step                           %f\n", dt); 
    printf(" time filter parameter               %f\n", alpha); 
  }

  // Start timer
  tstart = wtime(); 
  time = 0.;
  t100 = 0.;
  t200 = 0.;
  t300 = 0.;

  // ** Start of time loop ** 
#pragma omp parallel default (shared) private(i,j,ncycle,tdts8,tdtsdx,tdtsdy) firstprivate(tdt)

  for (ncycle=1;ncycle<=ITMAX;ncycle++) {
    
    // Compute capital u, capital v, z and h
    #pragma omp master
    c1 = wtime();

#pragma omp for schedule (static,chunk_size) nowait
    for (i=0;i<M;i++) {
      for (j=0;j<N;j++) {
        cu[i + 1][j] = .5 * (p[i + 1][j] + p[i][j]) * u[i + 1][j];
      }
    }

#pragma omp for schedule (static,chunk_size) nowait
    for (i=0;i<M;i++) {
      for (j=0;j<N;j++) {
        cv[i][j + 1] = .5 * (p[i][j + 1] + p[i][j]) * v[i][j + 1];
      }
    }

#pragma omp for schedule (static,chunk_size) nowait
    for (i=0;i<M;i++) {
      for (j=0;j<N;j++) {
        z[i + 1][j + 1] = (fsdx * (v[i + 1][j + 1] - v[i][j + 1]) - fsdy * (u[i + 1][j + 1] - u[i + 1][j])) / (p[i][j] + p[i + 1][j] + p[i + 1][j + 1] + p[i][j + 1]);
      }
    }

#pragma omp for schedule (static,chunk_size)
    for (i=0;i<M;i++) {
      for (j=0;j<N;j++) {
        h[i][j] = p[i][j] + .25 * (u[i + 1][j] * u[i + 1][j] + u[i][j] * u[i][j] + v[i][j + 1] * v[i][j + 1] + v[i][j] * v[i][j]);
      }
    }

    #pragma omp master
    {
        c2 = wtime();
        t100 = t100 + (c2 - c1); 
    }

    // Periodic continuation
    #pragma omp single
    {
        for (j=0;j<N;j++) {
          cu[0][j] = cu[M][j];
          cv[M][j + 1] = cv[0][j + 1];
          z[0][j + 1] = z[M][j + 1];
          h[M][j] = h[0][j];
        }
        cu[0][N] = cu[M][0];
        cv[M][0] = cv[0][N];
        z[0][0] = z[M][N];
        h[M][N] = h[0][0];
    }
#pragma omp for schedule (static,chunk_size)
    for (i=0;i<M;i++) {
      cu[i + 1][N] = cu[i + 1][0];
      cv[i][0] = cv[i][N];
      z[i + 1][0] = z[i + 1][N];
      h[i][N] = h[i][0];
    }

    // Compute new values u,v and p
    tdts8 = tdt / 8.;
    tdtsdx = tdt / dx;
    tdtsdy = tdt / dy;

    #pragma omp master
    c1 = wtime(); 

#pragma omp for schedule (static,chunk_size) nowait
    for (i=0;i<M;i++) {
      for (j=0;j<N;j++) {
        unew[i + 1][j] = uold[i + 1][j] + tdts8 * (z[i + 1][j + 1] + z[i + 1][j]) * (cv[i + 1][j + 1] + cv[i][j + 1] + cv[i][j] + cv[i + 1][j]) - tdtsdx * (h[i + 1][j] - h[i][j]);
      }
    }

#pragma omp for schedule (static,chunk_size) nowait
    for (i=0;i<M;i++) {
      for (j=0;j<N;j++) {
        vnew[i][j + 1] = vold[i][j + 1] - tdts8 * (z[i + 1][j + 1] + z[i][j + 1]) * (cu[i + 1][j + 1] + cu[i][j + 1] + cu[i][j] + cu[i + 1][j]) - tdtsdy * (h[i][j + 1] - h[i][j]);
      }
    }

#pragma omp for schedule (static,chunk_size)
    for (i=0;i<M;i++) {
      for (j=0;j<N;j++) {
        pnew[i][j] = pold[i][j] - tdtsdx * (cu[i + 1][j] - cu[i][j]) - tdtsdy * (cv[i][j + 1] - cv[i][j]); 
      }
    }

    #pragma omp master
    {
        c2 = wtime();  
        t200 = t200 + (c2 - c1); 
    }

    // Periodic continuation
    #pragma omp single
    {
        for (j=0;j<N;j++) {
          unew[0][j] = unew[M][j];
          vnew[M][j + 1] = vnew[0][j + 1];
          pnew[M][j] = pnew[0][j];
        }

        unew[0][N] = unew[M][0];
        vnew[M][0] = vnew[0][N];
        pnew[M][N] = pnew[0][0];
    }
#pragma omp for schedule (static,chunk_size)
    for (i=0;i<M;i++) {
      unew[i + 1][N] = unew[i + 1][0];
      vnew[i][0] = vnew[i][N];
      pnew[i][N] = pnew[i][0];
    }
    #pragma omp master
    time = time + dt;

    // Time smoothing and update for next cycle
    if ( ncycle > 1 ) {

      #pragma omp master
      c1 = wtime(); 

#pragma omp for schedule (static,chunk_size) nowait
      for (i=0;i<M_LEN;i++) {
        for (j=0;j<N_LEN;j++) {
          uold[i][j] = u[i][j] + alpha * (unew[i][j] - 2. * u[i][j] + uold[i][j]);
        }
      }

#pragma omp for schedule (static,chunk_size) nowait
      for (i=0;i<M_LEN;i++) {
        for (j=0;j<N_LEN;j++) {
          vold[i][j] = v[i][j] + alpha * (vnew[i][j] - 2. * v[i][j] + vold[i][j]);
        }
      }

#pragma omp for schedule (static,chunk_size) nowait
      for (i=0;i<M_LEN;i++) {
        for (j=0;j<N_LEN;j++) {
          pold[i][j] = p[i][j] + alpha * (pnew[i][j] - 2. * p[i][j] + pold[i][j]);
        }
      }

#pragma omp for schedule (static,chunk_size) nowait
      for (i=0;i<M_LEN;i++) {
        for (j=0;j<N_LEN;j++) {
          u[i][j] = unew[i][j];
        }
      }

#pragma omp for schedule (static,chunk_size) nowait
      for (i=0;i<M_LEN;i++) {
        for (j=0;j<N_LEN;j++) {
          v[i][j] = vnew[i][j];
        }
      }

#pragma omp for schedule (static,chunk_size) 
      for (i=0;i<M_LEN;i++) {
        for (j=0;j<N_LEN;j++) {
          p[i][j] = pnew[i][j];
        }
      }

      #pragma omp master
      {
          c2 = wtime(); 
          t300 = t300 + (c2 - c1);
      } 
    } else {

      tdt = tdt + tdt;

#pragma omp for schedule (static,chunk_size)
      for (i=0;i<M_LEN;i++) {
        for (j=0;j<N_LEN;j++) {
          uold[i][j] = u[i][j];
          vold[i][j] = v[i][j];
          pold[i][j] = p[i][j];
          u[i][j] = unew[i][j];
          v[i][j] = vnew[i][j];
          p[i][j] = pnew[i][j];
        }
      }
    }
  }

  // ** End of time loop ** 

  fprintf(stdout, "P CHECKSUM after %d steps = %15.7e\n", 
	  ITMAX, compute_checksum(pnew,M_LEN,N_LEN));
  fprintf(stdout, "U CHECKSUM after %d steps = %15.7e\n", 
	  ITMAX, compute_checksum(unew,M_LEN,N_LEN));
  fprintf(stdout, "V CHECKSUM after %d steps = %15.7e\n", 
	  ITMAX, compute_checksum(vnew,M_LEN,N_LEN));

  // Output p, u, v fields and run times.
  if (L_OUT) {
    c2 = wtime(); 
    ctime = c2 - tstart;
    tcyc = ctime / ITMAX;

    fprintf(stdout,"\n");
    fprintf(stdout," Job run on %d threads with a chunk size of %d\n",
	    nthreads, chunk_size);
    fprintf(stdout," No. of steps = %d, total time = %f, time per cycle = %f (s)\n", 
	    ITMAX, ctime, tcyc);
    fprintf(stdout," time for c{u,v},z,h calc = %.6f s\n", t100);
    fprintf(stdout," time for {u,v,p}new calc = %.6f s\n", t200);
    fprintf(stdout," time for time-smoothing  = %.6f s\n", t300);
  }

  return(0);
}
