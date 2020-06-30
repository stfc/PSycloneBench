#ifndef _KERNEL_BUILDER_INCLUDE
#define _KERNEL_BUILDER_INCLUDE

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#define STD_STRING_LEN 128

void init_device(
         int platform_selection,
         cl_device_id *device,
		 char *version_str,
		 cl_context *context);

const char* OCL_GetErrorString(cl_int error);

void check_status(const char *text, cl_int err);

cl_program get_program(cl_context context,
		       const cl_device_id *device,
		       int is_source_file,
		       const char *filename);

cl_program get_source_kernel(cl_context context,
			     const cl_device_id *device,
			     const char *filename);

cl_program get_binary_kernel(cl_context context,
			     const cl_device_id *device,
			     const char *filename);

cl_ulong duration_ns(cl_event event);
#endif
