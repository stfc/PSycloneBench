/** Contains prototypes for the C versions of the Momentum kernels */
#ifndef _MOM_KERNEL_INCLUDE
#define _MOM_KERNEL_INCLUDE

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

void set_args_momu(cl_kernel kern,
		   cl_int *nx,
		   cl_mem *ua_device,
		   cl_mem *un_device,
		   cl_mem *vn_device,
		   cl_mem *hu_device,
		   cl_mem *hv_device,
		   cl_mem *ht_device,
		   cl_mem *ssha_u_device,
		   cl_mem *sshn_device,
		   cl_mem *sshn_u_device,
		   cl_mem *sshn_v_device,
		   cl_mem *tmask_device,
		   cl_mem *e1u_device, cl_mem *e1v_device,
		   cl_mem *e1t_device, cl_mem *e2u_device,
		   cl_mem *e2t_device, cl_mem *e12u_device,
		   cl_mem *gphiu_device,
		   cl_double *rdt, cl_double *cbfr, cl_double *visc);

void set_args_momv(cl_kernel kern,
		   cl_int *nx,
		   cl_mem *va_device,
		   cl_mem *un_device,
		   cl_mem *vn_device,
		   cl_mem *hu_device,
		   cl_mem *hv_device,
		   cl_mem *ht_device,
		   cl_mem *ssha_v_device,
		   cl_mem *sshn_device,
		   cl_mem *sshn_u_device,
		   cl_mem *sshn_v_device,
		   cl_mem *tmask_device,
		   cl_mem *e1v_device, cl_mem *e1t_device,
		   cl_mem *e2u_device, cl_mem *e2v_device,
		   cl_mem *e2t_device, cl_mem *e12v_device,
		   cl_mem *gphiv_device,
		   cl_double *rdt, cl_double *cbfr, cl_double *visc);

void momentum_u_code(int ji, int jj, int width,
		     double *ua, double *un, double *vn,
		     double *hu, double *hv, double *ht, double *ssha_u,
		     double *sshn, double *sshn_u, double *sshn_v,
		     int *tmask,
		     double *e1u, double *e1v, double *e1t,
		     double *e2u, double *e2t, double *e12u, double *gphiu,
		     double rdt, double cbfr, double visc);

void momentum_v_code(int ji, int jj, int width,
		     double *va, double *un, double *vn, 
		     double *hu, double *hv, double *ht, double *ssha_v, 
		     double *sshn, double *sshn_u, double *sshn_v, 
		     int *tmask, double *e1v, double *e1t, double *e2u,
		     double *e2v, double *e2t, double *e12v, double *gphiv,
		     double rdt, double cbfr, double visc);
#endif
