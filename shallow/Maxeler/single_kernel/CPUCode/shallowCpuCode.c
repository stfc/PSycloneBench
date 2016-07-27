#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

double compute_checksum(int minx, int maxx, int miny, int maxy, double** field) {
  int i, j;
  double sum = 0.0;

  for (i = minx; i < maxx; ++i) {
    for (j = miny; j < maxy; ++j) {
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
      printf("%e ", array[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void free_array(double** array) {
  free(array[0]);
  free(array);
}

int main(void)
{
  int m = 15, n = 15; // Global domain size
  int itmax = 1; // Number of timesteps
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

    printf("psi initial CHECKSUM = %e\n", compute_checksum(1, m_len, 1, n_len, psi));
    printf("p initial CHECKSUM = %e\n", compute_checksum(0, m, 0, n, p));
    printf("u initial CHECKSUM = %e\n", compute_checksum(1, m_len, 0, n, u));
    printf("v initial CHECKSUM = %e\n", compute_checksum(0, m, 1, n_len, v));
  }

  max_file_t *maxfile = shallow_init();
  max_engine_t *engine = max_load(maxfile, "*");

  shallow_actions_t run_scalar;
  run_scalar.param_len = m_len * n_len;
  run_scalar.param_fsdx = fsdx;
  run_scalar.param_fsdy = fsdy;
  run_scalar.param_tdt = tdt;
  run_scalar.param_dx = dx;
  run_scalar.param_dy = dy;
  run_scalar.instream_p = p[0];
  run_scalar.instream_u = u[0];
  run_scalar.instream_v = v[0];
  run_scalar.instream_pold = pold[0];
  run_scalar.instream_uold = uold[0];
  run_scalar.instream_vold = vold[0];
  run_scalar.outstream_unew = unew[0];
  run_scalar.outstream_vnew = vnew[0];
  run_scalar.outstream_pnew = pnew[0];

  tstart = wtime();

  shallow_run(engine, &run_scalar);

   c2 = wtime();
   ctime = c2 - tstart;
   tcyc = ctime / itmax;

   printf("\n");
   printf(" No. of steps = %d, total time = %f, time per cycle = %f (s)\n",
           itmax, ctime, tcyc);

  tdt = tdt + tdt;
  for (i = 0; i < m_len; i++) {
    for (j = 0; j < n_len; j++) {
      u[i][j] = unew[i][j];
      v[i][j] = vnew[i][j];
      p[i][j] = pnew[i][j];
    }
   }

  // End of time loop
  printf("P CHECKSUM after %d steps = %e\n",
          itmax, compute_checksum(0, m, 0, n, pnew));
  printf("U CHECKSUM after %d steps = %e\n",
          itmax, compute_checksum(1, m_len, 0, n, unew));
  printf("V CHECKSUM after %d steps = %e\n",
          itmax, compute_checksum(0, m, 1, n_len, vnew));


  max_unload(engine);

  printf("Done.\n");
  return 0;
}
