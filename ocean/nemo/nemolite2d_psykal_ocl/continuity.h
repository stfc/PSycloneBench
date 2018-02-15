#ifndef __CONTINUITY_HEADER
#define __CONTINUITY_HEADER

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

void set_args_continuity(cl_kernel cont_kernel,
			 cl_int *nx,
			 cl_mem *ssha_device,
			 cl_mem *sshn_device,
			 cl_mem *sshn_u_device,
			 cl_mem *sshn_v_device,
			 cl_mem *hu_device,
			 cl_mem *hv_device,
			 cl_mem *un_device,
			 cl_mem *vn_device,
			 cl_double *rdt,
			 cl_mem *e12t_device);

void continuity_code(int ji, int jj, int width,                     
		     double *ssha,
		     double *sshn,
		     double *sshn_u,
		     double *sshn_v,
		     double* hu,
		     double *hv,
		     double *un,
		     double *vn,
		     double rdt,
		     double *e12t);
#endif
