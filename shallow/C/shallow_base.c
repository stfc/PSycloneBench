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
//! xlf90 -qmaxmem=-1 -g -o shallow4 -qfixed=132 -qsclk=micro
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
/* Simplified C version of shallow based on both the OpenMP C and
 * simplified Fortran versions
 * 23 June 2016
 * */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

double compute_checksum(int minx, int maxx, int miny, int maxy, double** field) {
  int i, j;
  double sum = 0.0;

  for (j = miny; j < maxy; ++j) {
    for (i = minx; i < maxx; ++i) {
      sum += fabs(field[i][j]);
    }
  }

  return sum;
}

double wtime() {
   static int sec = -1;
   struct timeval tv;
   gettimeofday(&tv, NULL);
   if (sec < 0) sec = tv.tv_sec;
   return (tv.tv_sec - sec) + 1.0e-6 * tv.tv_usec;
}

double** allocate_array(int rows, int cols) {
  double **array = malloc(rows * sizeof(double*));
  array[0] = malloc(rows * cols * sizeof(double));
  for (int i = 1; i < rows; ++i) {
    array[i] = array[0] + i * cols;
  }
  return array;
}

void print_array(int rows, int cols, double** array) {
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      printf("%.16g", array[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void free_array(double** array) {
  free(array[0]);
  free(array);
}

int main() {
  int m = 128, n = 128; // Global domain size
  int itmax = 2000; // Number of timesteps
  bool l_out = true; // Produce output
  int m_len = m + 1, n_len = n + 1;

  double **u, **v, **p;
  double **unew, **vnew, **pnew;
  double **uold, **vold, **pold;
  double **cu, **cv, **z, **h, **psi;

  // Allocate the arrays
  u = allocate_array(m_len, n_len);
  v = allocate_array(m_len, n_len);
  p = allocate_array(m_len, n_len);
  unew = allocate_array(m_len, n_len);
  vnew = allocate_array(m_len, n_len);
  pnew = allocate_array(m_len, n_len);
  uold = allocate_array(m_len, n_len);
  vold = allocate_array(m_len, n_len);
  pold = allocate_array(m_len, n_len);
  cu = allocate_array(m_len, n_len);
  cv = allocate_array(m_len, n_len);
  z = allocate_array(m_len, n_len);
  h = allocate_array(m_len, n_len);
  psi = allocate_array(m_len, n_len);

  double dt, tdt, dx, dy, a, alpha, el, pi, tpi, di, dj, pcf;
  double tdts8, tdtsdx, tdtsdy, fsdx, fsdy;

  int ncycle;
  int i, j;

  // timer variables
  double t100, t200, t300;
  double tstart, ctime, tcyc;
  double c1, c2;

  dt = 90;
  tdt = dt; // Only for the first cycle

  dx = 100000;
  dy = 100000;
  fsdx = 4.0 / dx;
  fsdy = 4.0 / dy;


  a = 1000000;
  alpha = 0.001;

  el = n * dx;
  pi = 4 * atan(1.0);
  tpi = pi + pi;
  di = tpi / m;
  dj = tpi / n;
  pcf = pi * pi * a * a / (el * el);

  // Initial values of the stream functions and p
  for (i = 0; i < m_len; i++) {
    for (j = 0; j < n_len; j++) {
      psi[i][j] = a * sin((i + 0.5) * di) * sin((j + 0.5) * dj);
      p[i][j] = pcf * (cos(2.0 * (i) * di) + cos(2.0 * (j) * dj)) + 50000.0;
    }
  }

  // Initialize velocities
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
      u[i + 1][j] = -(psi[i + 1][j + 1] - psi[i + 1][j]) / dy;
      v[i][j + 1] = (psi[i + 1][j + 1] - psi[i][j + 1]) / dx;
    }
  }

  // Periodic continuation
  for (j = 0; j < n; ++j) {
    u[0][j] = u[m][j];
    v[m][j + 1] = v[0][j + 1];
  }
  for (i = 0; i < m; ++i) {
    u[i + 1][n] = u[i + 1][0];
    v[i][0] = v[i][n];
  }
  u[0][n] = u[m][0];
  v[m][0] = v[0][n];

  for (i=0; i < m_len; i++) {
    for (j = 0; j < n_len; j++) {
      uold[i][j] = u[i][j];
      vold[i][j] = v[i][j];
      pold[i][j] = p[i][j];
    }
  }

  // Print initial values
  if ( l_out ) {
    printf(" number of points in the x direction %d\n", n);
    printf(" number of points in the y direction %d\n", m);
    printf(" grid spacing in the x direction     %f\n", dx);
    printf(" grid spacing in the y direction     %f\n", dy);
    printf(" time step                           %f\n", dt);
    printf(" time filter parameter               %f\n", alpha);

    printf("psi initial CHECKSUM = %.16g\n", compute_checksum(1, m_len, 1, n_len, psi));
    printf("p initial CHECKSUM = %.16g\n", compute_checksum(0, m, 0, n, p));
    printf("u initial CHECKSUM = %.16gn", compute_checksum(1, m_len, 0, n, u));
    printf("v initial CHECKSUM = %.16g\n", compute_checksum(0, m, 1, n_len, v));
  }

  tstart = wtime();
  t100 = t200 = t300 = 0;

  for (ncycle = 1; ncycle <= itmax; ++ncycle) {
    c1 = wtime();

    // Compute capital u, capital v, z and h
    for (i = 0; i < m; ++i) {
       for (j = 0; j < n; ++j) {
         cu[i + 1][j] = 0.5 * (p[i + 1][j] + p[i][j]) * u[i + 1][j];
         cv[i][j + 1] = 0.5 * (p[i][j + 1] + p[i][j]) * v[i][j + 1];
         z[i + 1][j + 1] = (fsdx * (v[i + 1][j + 1] - v[i][j + 1]) - fsdy *
                           (u[i + 1][j + 1] - u[i + 1][j])) / (p[i][j] +
                            p[i + 1][j] + p[i + 1][j + 1] + p[i][j + 1]);
         h[i][j] = p[i][j] + 0.25 * (u[i + 1][j] * u[i + 1][j] + u[i][j] *
                   u[i][j] + v[i][j + 1] * v[i][j + 1] + v[i][j] * v[i][j]);
      }
    }

    c2 = wtime();
    t100 = t100 + (c2 - c1);

    // Periodic continuation
    for (j = 0; j < n; ++j) {
    	cu[0][j] = cu[m][j];
        cv[m][j + 1] = cv[0][j + 1];
        z[0][j + 1] = z[m][j + 1];
        h[m][j] = h[0][j];
    }
    for (i = 0; i < m; ++i) {
      cu[i + 1][n] = cu[i + 1][0];
      cv[i][0] = cv[i][n];
      z[i + 1][0] = z[i + 1][n];
      h[i][n] = h[i][0];
    }
    cu[0][n] = cu[m][0];
    cv[m][0] = cv[0][n];
    z[0][0] = z[m][n];
    h[m][n] = h[0][0];

    // Compute new values u,v and p
    tdts8 = tdt / 8.0;
    tdtsdx = tdt / dx;
    tdtsdy = tdt / dy;

    c1 = wtime();

    for (i = 0; i < m; ++i) {
      for (j = 0; j < n; ++j) {
        unew[i + 1][j] = uold[i + 1][j] + tdts8 * (z[i + 1][j + 1] + z[i + 1][j]) *
                         (cv[i + 1][j + 1] + cv[i][j + 1] + cv[i][j] + cv[i + 1][j]) -
                         tdtsdx * (h[i + 1][j] - h[i][j]);
        vnew[i][j + 1] = vold[i][j + 1] - tdts8 * (z[i + 1][j + 1] + z[i][j + 1]) *
                         (cu[i + 1][j + 1] + cu[i][j + 1] + cu[i][j] + cu[i + 1][j]) -
                         tdtsdy * (h[i][j + 1] - h[i][j]);
        pnew[i][j] = pold[i][j] - tdtsdx * (cu[i + 1][j] - cu[i][j]) - tdtsdy *
                     (cv[i][j + 1] - cv[i][j]);
      }
    }

    c2 = wtime();
    t200 = t200 + (c2 - c1);

    // Periodic continuation
    for (j = 0; j < n; ++j) {
      unew[0][j] = unew[m][j];
      vnew[m][j + 1] = vnew[0][j + 1];
      pnew[m][j] = pnew[0][j];
    }
    for (i = 0; i < m; ++i) {
      unew[i + 1][n] = unew[i + 1][0];
      vnew[i][0] = vnew[i][n];
      pnew[i][n] = pnew[i][0];
    }
    unew[0][n] = unew[m][0];
    vnew[m][0] = vnew[0][n];
    pnew[m][n] = pnew[0][0];

    // Time smoothing and update for next cycle
    if (ncycle > 1) {
      c1 = wtime();
      for (i = 0; i < m_len; ++i) {
        for (j = 0; j < n_len; ++j) {
          uold[i][j] = u[i][j] + alpha * (unew[i][j] - 2.0 * u[i][j] + uold[i][j]);
          vold[i][j] = v[i][j] + alpha * (vnew[i][j] - 2.0 * v[i][j] + vold[i][j]);
          pold[i][j] = p[i][j] + alpha * (pnew[i][j] - 2.0 * p[i][j] + pold[i][j]);
          u[i][j] = unew[i][j];
          v[i][j] = vnew[i][j];
          p[i][j] = pnew[i][j];
        }
      }
      c2 = wtime();
      t300 = t300 + (c2 - c1);
    } else {
      tdt = tdt + tdt;
      for (i = 0; i < m_len; i++) {
        for (j = 0; j < n_len; j++) {
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

  // End of time loop
  printf("P CHECKSUM after %d steps = %.16g\n",
          itmax, compute_checksum(0, m, 0, n, pnew));
  printf("U CHECKSUM after %d steps = %.16g\n",
          itmax, compute_checksum(1, m_len, 0, n, unew));
  printf("V CHECKSUM after %d steps = %.16g\n",
          itmax, compute_checksum(0, m, 1, n_len, vnew));

  // Output p, u, v fields and run times.
  if (l_out) {
    c2 = wtime();
    ctime = c2 - tstart;
    tcyc = ctime / itmax;

    printf("\n");
    printf(" No. of steps = %d, total time = %f, time per cycle = %f (s)\n",
            itmax, ctime, tcyc);
    printf(" time for c{u,v},z,h calc = %.6f s\n", t100);
    printf(" time for {u,v,p}new calc = %.6f s\n", t200);
    printf(" time for time-smoothing  = %.6f s\n", t300);
  }

  // Free the arrays
  free_array(u);
  free_array(v);
  free_array(p);
  free_array(unew);
  free_array(vnew);
  free_array(pnew);
  free_array(uold);
  free_array(vold);
  free_array(pold);
  free_array(cu);
  free_array(cv);
  free_array(z);
  free_array(h);
  free_array(psi);

  return 0;
}
