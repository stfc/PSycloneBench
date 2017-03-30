#include "continuity.h"
#include <stdlib.h>
#include <stdio.h>

int main(){
  int nx = 128;
  int ny = 128;
  int ji, jj, idx;
  double *ssha, *sshn, *sshn_u, *sshn_v;
  double *hu, *hv, *un, *vn;
  double rdt;
  double *e12t;
  double dep_const = 2.0;
  double dx = 0.5, dy = 0.5;
  ji = 1;
  jj = 1;
  rdt = 0.5;

  /* Field initialisation */
  
  ssha = malloc(nx*ny*sizeof(double));
  sshn = malloc(nx*ny*sizeof(double));
  sshn_u = malloc(nx*ny*sizeof(double));
  sshn_v = malloc(nx*ny*sizeof(double));
  hu = malloc(nx*ny*sizeof(double));
  hv = malloc(nx*ny*sizeof(double));
  un = malloc(nx*ny*sizeof(double));
  vn = malloc(nx*ny*sizeof(double));
  e12t = malloc(nx*ny*sizeof(double));

  for(jj=0;jj<ny;jj++){
    for(ji=0;ji<nx;ji++){
      idx = jj*nx + ji;
      hu[idx] = dep_const;
      hv[idx] = dep_const;
      un[idx] = ji;
      vn[idx] = jj;
      sshn_u[idx] = 0.0;
      sshn_v[idx] = 0.0;
      sshn[idx] = 1.0;
      e12t[idx] = dx*dy;
      ssha[idx] = 0.0;
    }
  }

  /* Run the kernel */
  for(jj=1;jj<ny;jj++){
    for(ji=1;ji<nx;ji++){
      continuity_code(ji,jj, nx,                     
		      ssha, sshn, sshn_u, sshn_v, hu, hv,
		      un, vn, rdt, e12t);
    }
  }

  /* Compute and output a checksum */
  double *sshaptr = &(ssha[0]);
  double sum = 0.0;
  for(jj=0;jj<ny;jj++){
    for(ji=0;ji<nx;ji++){
      sum += *(sshaptr++);
    }
  }
  fprintf(stdout, "Checksum = %e\n", &sum);
}

