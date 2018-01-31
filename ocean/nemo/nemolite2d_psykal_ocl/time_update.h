#ifndef _TIME_UPDATE_INCLUDE
#define _TIME_UPDATE_INCLUDE

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

void next_sshu_code(int ji, int jj, int width,
		    double* sshn_u,
		    double* sshn,
		    int* tmask,
		    double* e12t,
		    double* e12u);

void set_args_next_sshu(cl_kernel kern,
			cl_int *nx,
			cl_mem *sshn_u_device,
			cl_mem *sshn_device,
			cl_mem *tmask_device,
			cl_mem *e12t_device,
			cl_mem *e12u_device);

void set_args_next_sshv(cl_kernel kern,
			cl_int *nx,
			cl_mem* sshn_v_device,
			cl_mem* sshn_device,
			cl_mem* tmask_device,
			cl_mem* e12t_device,
			cl_mem* e12v_device);

void next_sshv_code(int ji, int jj, int width,
		    double* sshn_v,
		    double* sshn,
		    int* tmask,
		    double* e12t,
		    double* e12v);

#endif
