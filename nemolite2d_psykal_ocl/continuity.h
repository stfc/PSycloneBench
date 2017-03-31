#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

__kernel void continuity_code(int width,                     
		     __global double *ssha,
		     __global double *sshn,
		     __global double *sshn_u,
		     __global double *sshn_v,
		     __global double* hu,
		     __global double *hv,
		     __global double *un,
		     __global double *vn,
		     double rdt,
		     __global double *e12t);
