// OpenCL set_args definition
#ifdef OPENCL_HOST

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

void set_args_continuity(cl_kernel cont_kernel,
			 cl_int *nx,
			 cl_int *xstart,
			 cl_int *xstop,
			 cl_int *ystart,
			 cl_int *ystop,
			 cl_mem *ssha_device,
			 cl_mem *sshn_device,
			 cl_mem *sshn_u_device,
			 cl_mem *sshn_v_device,
			 cl_mem *hu_device,
			 cl_mem *hv_device,
			 cl_mem *un_device,
			 cl_mem *vn_device,
			 cl_double *rdt,
			 cl_mem *e12t_device){
  cl_int ret;
  int arg_idx = 0;
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_int), (void *)nx);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_int), (void *)xstart);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_int), (void *)xstop);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_int), (void *)ystart);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_int), (void *)ystop);
  check_status("clSetKernelArg", ret);

  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)ssha_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_u_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)sshn_v_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)hu_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)hv_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)un_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)vn_device);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_double),
		       (void *)rdt);
  check_status("clSetKernelArg", ret);
  ret = clSetKernelArg(cont_kernel, arg_idx++, sizeof(cl_mem),
		       (void *)e12t_device);
  check_status("clSetKernelArg", ret);
  
  fprintf(stdout, "Set %d arguments for Continuity kernel\n", arg_idx);
}
#endif

// OpenCL Kernel definition
#ifdef __OPENCL_VERSION__
__kernel void continuity_code(
        int width,
        int xstart, int xstop,
        int ystart, int ystop,
        __global double* restrict ssha,
        __global double* restrict sshn,
        __global double* restrict sshn_u,
        __global double* restrict sshn_v,
        __global double* restrict hu,
        __global double* restrict hv,
        __global double* restrict un,
        __global double* restrict vn,
        double rdt,
        __global double* restrict e12t
){
    for(int jj = ystart; jj <= ystop; jj++){
        for(int ji = xstart; ji <= xstop; ji++){
            int idx = jj * width + ji;
            int idxim1 = idx - 1;
            int idxjm1 = idx - width;

            const double rtmp1 = (sshn_u[idx] + hu[idx]) * un[idx];
            const double rtmp2 = (sshn_u[idxim1] + hu[idxim1]) * un[idxim1];
            const double rtmp3 = (sshn_v[idx] + hv[idx]) * vn[idx];
            const double rtmp4 = (sshn_v[idxjm1] + hv[idxjm1]) * vn[idxjm1];

            ssha[idx] = sshn[idx] + (rtmp2 - rtmp1 + rtmp4 - rtmp3) *
              rdt / e12t[idx];
        }
    }
}
#endif
