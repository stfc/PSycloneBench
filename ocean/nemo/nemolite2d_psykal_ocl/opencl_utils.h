#ifndef _KERNEL_BUILDER_INCLUDE
#define _KERNEL_BUILDER_INCLUDE

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#define STD_STRING_LEN 128

void init_device(cl_device_id *device,
		 char *version_str,
		 cl_context *context);

const char* OCL_GetErrorString(cl_int error);

void check_status(const char *text, cl_int err);

cl_kernel get_kernel(cl_context *context,
		     cl_device_id *device,
		     const char *version_str,
		     const char *filename,
		     const char *kernel_name);

cl_program get_source_kernel(cl_context *context,
			     cl_device_id *device,
			     const char *filename);

cl_program get_binary_kernel(cl_context *context,
			     cl_device_id *device,
			     const char *filename);

cl_ulong duration_ns(cl_event event);
#endif
