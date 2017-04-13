#ifndef _KERNEL_BUILDER_INCLUDE
#define _KERNEL_BUILDER_INCLUDE

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

const char* OCL_GetErrorString(cl_int error);

void check_status(char *text, cl_int err);

cl_kernel build_kernel(cl_context *context, cl_device_id *device,
		       char *filename, char *kernel_name);

#endif
