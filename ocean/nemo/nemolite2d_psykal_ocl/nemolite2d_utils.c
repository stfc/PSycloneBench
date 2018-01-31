#include "opencl_utils.h"
#include "math.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

int init_fields(int nx, int ny, cl_double dx, cl_double dy,
		cl_double dep_const,
		int *xstart, int *xstop, int *ystart, int *ystop,
		cl_double **ssha, cl_double **ssha_u, cl_double **ssha_v,
		cl_double **sshn, cl_double **sshn_u, cl_double **sshn_v,
		cl_double **hu, cl_double **hv, cl_double **ht,
		cl_double **un, cl_double **vn, cl_double **ua, cl_double **va,
		cl_double **gphiu, cl_double **gphiv,
		cl_double **e1u, cl_double **e1v, cl_double **e1t,
		cl_double **e2u, cl_double **e2v, cl_double **e2t,
		cl_double **e12t, cl_double **e12u, cl_double **e12v,
		cl_int **tmask){
  int ji, jj;
  int idx;
  long buff_size = nx * ny * sizeof(cl_double);
  
  *xstart = 1;
  *xstop = nx - 1;
  *ystart = 1;
  *ystop = ny - 1;

  *ssha = (cl_double*)malloc(buff_size);
  *ssha_u = (cl_double*)malloc(buff_size);
  *ssha_v = (cl_double*)malloc(buff_size);
  *sshn = (cl_double*)malloc(buff_size);
  *sshn_u = (cl_double*)malloc(buff_size);
  *sshn_v = (cl_double*)malloc(buff_size);
  *hu = (cl_double*)malloc(buff_size);
  *hv = (cl_double*)malloc(buff_size);
  *ht = (cl_double*)malloc(buff_size);
  *un = (cl_double*)malloc(buff_size);
  *vn = (cl_double*)malloc(buff_size);
  *ua = (cl_double*)malloc(buff_size);
  *va = (cl_double*)malloc(buff_size);
  *e1u = (cl_double*)malloc(buff_size);
  *e1v = (cl_double*)malloc(buff_size);
  *e1t = (cl_double*)malloc(buff_size);
  *e2u = (cl_double*)malloc(buff_size);
  *e2v = (cl_double*)malloc(buff_size);
  *e2t = (cl_double*)malloc(buff_size);
  *e12u = (cl_double*)malloc(buff_size);
  *e12v = (cl_double*)malloc(buff_size);
  *e12t = (cl_double*)malloc(buff_size);
  *tmask = (cl_int*)malloc(nx*ny*sizeof(cl_int));

  *gphiu = (cl_double*)malloc(buff_size);
  *gphiv = (cl_double*)malloc(buff_size);

  for(jj=0;jj<ny;jj++){
    for(ji=0;ji<nx;ji++){
      idx = jj*nx + ji;
      (*hu)[idx] = dep_const;
      (*hv)[idx] = dep_const;
      (*ht)[idx] = dep_const;
      (*un)[idx] = 0.01;
      (*vn)[idx] = 0.0;
      (*sshn_u)[idx] = cos((360.0*ji)/(float)nx);
      (*sshn_v)[idx] = cos((360.0*ji)/(float)nx);
      (*sshn)[idx] = sin((360.0*ji)/(float)nx);
      (*ssha)[idx] = 0.0;
      // Grid properties
      (*e1u)[idx] = dx;
      (*e1v)[idx] = dx;
      (*e1t)[idx] = dx;
      (*e2u)[idx] = dy;
      (*e2v)[idx] = dy;
      (*e2t)[idx] = dy;
      (*e12u)[idx] = dx*dy;
      (*e12v)[idx] = (*e12u)[idx];
      (*e12t)[idx] = (*e12u)[idx];
      // f-plane test case (constant Coriolis parameter)
      (*gphiu)[idx] = 50.0;
      (*gphiv)[idx] = 50.0;
      // All inner cells
      (*tmask)[idx] = 1;
    }
  }
  for(jj=0;jj<ny;jj++){
    idx = jj*nx;
    // West solid boundary
    for(ji=0; ji<*xstart; ji++){
      (*tmask)[idx+ji] = 0;
    }
    // East solid boundary
    for(ji=*xstop; ji<nx; ji++){
      (*tmask)[idx+ji] = 0;
    }
  }
  // Southern open boundary
  for(jj=0; jj<*ystart; jj++){
    idx = jj*nx;
    for(ji=0;ji<nx;ji++){
      (*tmask)[idx + ji] = -1;
    }
  }
  // North solid boundary
  for(jj=*ystop; jj<ny; jj++){
    idx = jj*nx;
    for(ji=0;ji<nx;ji++){
      (*tmask)[idx + ji] = 0;
    }
  }

  return 0;
}
