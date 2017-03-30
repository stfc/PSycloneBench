#include "continuity.h"
#include <stdlib.h>
#include <stdio.h>

int main(){
  int nx = 128;
  int ny = 128;
  int ji, jj;
  double *ssha, *sshn, *sshn_u, *sshn_v;
  double *hu, *hv, *un, *vn;
  double rdt;
  double *e12t;
  double dep_const = 2.0;
  
  ji = 1;
  jj = 1;
  rdt = 0.5;

  ssha = malloc(nx*ny*sizeof(double));
  sshn = malloc(nx*ny*sizeof(double));
  sshn_u = malloc(nx*ny*sizeof(double));
  sshn_v = malloc(nx*ny*sizeof(double));
  hu = malloc(nx*ny*sizeof(double));
  hv = malloc(nx*ny*sizeof(double));
  un = malloc(nx*ny*sizeof(double));
  vn = malloc(nx*ny*sizeof(double));
  e12t = malloc(nx*ny*sizeof(double));

  double *huptr = &(hu[0]);
  double *hvptr = &(hv[0]);
  for(jj=0;jj<ny;jj++){
    for(ji=0;ji<nx;ji++){
      *(huptr++) = dep_const;
      *(hvptr++) = dep_const;
    }
  }
  
  for(jj=0;jj<ny;jj++){
    for(ji=0;ji<nx;ji++){
      continuity_code(ji,jj, nx,                     
		      ssha, sshn, sshn_u, sshn_v, hu, hv,
		      un, vn, rdt, e12t);
    }
  }

  double *sshaptr = &(ssha[0]);
  double sum = 0.0;
  for(jj=0;jj<ny;jj++){
    for(ji=0;ji<nx;ji++){
      sum += *(sshaptr++);
    }
  }
  fprintf(stdout, "Checksum = %e\n", &sum);
}

