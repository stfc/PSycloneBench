
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

int init_fields(int nx, int ny,  cl_double dx, cl_double dy,
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
		cl_int **tmask);
