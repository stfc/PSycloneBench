#include "opencl_utils.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

// Kernels
#define OPENCL_HOST
#include "../../kernels/c_family/continuity_kern.c"
#include "../../kernels/c_family/momentum_u_kern.c"
#include "../../kernels/c_family/momentum_v_kern.c"
#include "../../kernels/c_family/boundary_conditions_kern.c"
#include "../../kernels/c_family/time_update_kern.c"

/** Maximum number of OpenCL devices we will query */
#define MAX_DEVICES 4

/* Number of OpenCL command queues we will use */
#define NUM_QUEUES 2

int first_time = 1;

enum KERNELS {
  K_CONTINUITY,
  K_MOM_U,
  K_MOM_V,
  K_BC_SSH,
  K_BC_SOLID_U,
  K_BC_SOLID_V,
  K_BC_FLATHER_U,
  K_BC_FLATHER_V,
  K_NEXT_SSH_U,
  K_NEXT_SSH_V,
  K_NUM_KERNELS
};

/* Names of all of our kernels */
static const char* kernel_names[K_NUM_KERNELS] =
{
  "continuity_code",
  "momentum_u_code",
  "momentum_v_code",
  "bc_ssh_code",
  "bc_solid_u_code",
  "bc_solid_v_code",
  "bc_flather_u_code",
  "bc_flather_v_code",
  "next_sshu_code",
  "next_sshv_code"
};

/* The source files containing each kernel (only used if not compiling
   them to a single image) */
static const char* kernel_files[K_NUM_KERNELS] =
{
  "../../kernels/c_family/continuity_kern.c",
  "../../kernels/c_family/momentum_u_kern.c",
  "../../kernels/c_family/momentum_v_kern.c",
  "../../kernels/c_family/boundary_conditions_kern.c",
  "../../kernels/c_family/boundary_conditions_kern.c",
  "../../kernels/c_family/boundary_conditions_kern.c",
  "../../kernels/c_family/boundary_conditions_kern.c",
  "../../kernels/c_family/boundary_conditions_kern.c",
  "../../kernels/c_family/time_update_kern.c",
  "../../kernels/c_family/time_update_kern.c"   
};

// OpenCL variables should be permanent between multiple invoke calls
cl_device_id device;
cl_context context = NULL;
/* String holding info on chosen device */
char version_str[STD_STRING_LEN];

/* In order to run multiple kernels concurrently, we must have
 multiple command queues (since the Intel OpenCL SDK only supports
 in-order queues) */
cl_command_queue command_queue[NUM_QUEUES];
cl_command_queue_properties queue_properties = 0;
cl_bool profiling_enabled = 0;

/* Buffers on the device */
cl_mem ssha_device = NULL;
cl_mem ssha_u_device = NULL;
cl_mem ssha_v_device = NULL;
cl_mem sshn_device = NULL;
cl_mem sshn_u_device = NULL;
cl_mem sshn_v_device = NULL;
cl_mem ht_device = NULL;
cl_mem hu_device = NULL;
cl_mem hv_device = NULL;
cl_mem un_device = NULL;
cl_mem vn_device = NULL;
cl_mem ua_device = NULL;
cl_mem va_device = NULL;
cl_mem e1u_device = NULL;
cl_mem e1v_device = NULL;
cl_mem e1t_device = NULL;
cl_mem e2u_device = NULL;
cl_mem e2v_device = NULL;
cl_mem e2t_device = NULL;
cl_mem e12u_device = NULL;
cl_mem e12v_device = NULL;
cl_mem e12t_device = NULL;
cl_mem gphiu_device = NULL;
cl_mem gphiv_device = NULL;
cl_mem tmask_device = NULL;

cl_program program = NULL;

/* Our OpenCL kernel objects */
cl_kernel clkernel[K_NUM_KERNELS];
cl_event clkernevt[K_NUM_KERNELS];



void c_invoke_time_step(
        // Fields
        double * ssha_t,
        double * sshn_t,
        double * sshn_u,
        double * sshn_v,
        double * hu,
        double * hv,
        double * un,
        double * vn,
        double * ua,
        double * ht,
        double * ssha_u,
        double * va,
        double * ssha_v,
        // Grid
        int * tmask,
        double * area_t,
        double * area_u,
        double * area_v,
        double * dx_u,
        double * dx_v,
        double * dx_t,
        double * dy_u,
        double * dy_v,
        double * dy_t,
        double * gphiu,
        // Scalars
        int istep,
        int internal_xstart,
        int internal_xstop,
        int internal_ystart,
        int internal_ystop,
        int width,
        double rdt,
        double cbfr,
        double visc,
        double omega,
        double d2r,
        double g
        ){

    if(first_time){
        printf("OpenCL initialization\n");

        /* Run-time configuration of the benchmark */
        if(getenv("NEMOLITE2D_PROFILING")){
            profiling_enabled = CL_TRUE;
            /* We create the OpenCL command queue with profiling enabled */
            queue_properties = CL_QUEUE_PROFILING_ENABLE;
        }

        /* Check to see whether we should get our kernels from a single image file */
        char *env_string;
        char *image_file = NULL;
        if( (env_string = getenv("NEMOLITE2D_SINGLE_IMAGE")) ){
            /* the %ms below instructs sscanf to allocate memory as required */
            if(sscanf(env_string, "%ms", &image_file) != 1){
                fprintf(stderr, "Error parsing NEMOLITE2D_SINGLE_IMAGE environment "
                        "variable (%s)\n", env_string);
                exit(1);
            }
        }
        
        int platform = 0; // By default choose platform number 0
        if( (env_string = getenv("OPENCL_PLATFORM")) ){
            platform = atoi(env_string); // If not valid conversion it also returns 0
        }

        init_device(platform, &device, version_str, &context);

        int ret;
        for(int ji=0; ji<NUM_QUEUES; ji++){
            /* The Intel/Altera OpenCL SDK is only version 1.0 */
            /* NVIDIA only support OpenCL 1.2 so we get a seg. fault if we attempt
            to call the ...WithProperties version of this routine */
            command_queue[ji] = clCreateCommandQueue(
                    context, device, queue_properties, &ret);
            check_status("clCreateCommandQueue", ret);
        }

        // Create our Program object (contains all of the individual kernels)
        if(image_file){
            /* We expect all kernels to have been compiled into a single image */
            fprintf(stderr, " Single image version not supported at the moment.");
            exit(1);
            program = get_program(context, &device, version_str, image_file);
        }else{
            //fprintf(stderr, "Please set NEMOLITE2D_SINGLE_IMAGE to point to the "
            //                ".aocx file containing the compiled kernels\n");
            //exit(1);
        }
        
        /* Create OpenCL Kernels and associated event objects (latter used
        to obtain detailed timing information). */
        for(int ikern=0; ikern<K_NUM_KERNELS; ikern++){
            if(!image_file){
                program = get_program(context, &device, version_str,
                    kernel_files[ikern]);
            }
            fprintf(stdout, "Creating kernel %s...\n", kernel_names[ikern]);
            clkernel[ikern] = clCreateKernel(program, kernel_names[ikern], &ret);
            check_status("clCreateKernel", ret);
        } 
  
        /* Create Device Memory Buffers -- just works with square domains nx==ny */
        int num_buffers = 0;
        int buff_size = width*width*sizeof(cl_double);
        ssha_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			                         NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        ssha_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        ssha_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        sshn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
                                     NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        sshn_u_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
				                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        sshn_v_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
				                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        hu_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        hv_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        ht_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
  
        /* Velocity fields */
        un_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        vn_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
#ifndef SINGLE_KERNEL
        ua_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        va_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
#endif

        /* Mesh scale factors */
        e1t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        e1u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        e1v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        e2u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        e2v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        e2t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        e12u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                                     NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        e12v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                                     NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        tmask_device = clCreateBuffer(context, CL_MEM_READ_ONLY,
				        (size_t)(width*width*sizeof(cl_int)), NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        e12t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                                     NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);

        /* Coriolis parameters */
        gphiu_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
				                      NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        gphiv_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
				                      NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);

        fprintf(stdout, "Created %d device buffers OK\n", num_buffers);

        /* Set OpenCL Kernel Parameters for Continuity */
        set_args_continuity(clkernel[K_CONTINUITY],
		    &width,
		    &ssha_device, &sshn_device,
		    &sshn_u_device, &sshn_v_device,
		    &hu_device, &hv_device,
		    &un_device, &vn_device,
		    &rdt, &e12t_device); 
        /* Set OpenCL Kernel Parameters for Momentum-u */
        set_args_momu(clkernel[K_MOM_U],
            &width,
            &ua_device, &un_device, &vn_device,
            &hu_device, &hv_device, &ht_device,
            &ssha_u_device, &sshn_device,
            &sshn_u_device, &sshn_v_device,
            &tmask_device,
            &e1u_device, &e1v_device,
            &e1t_device, &e2u_device,
            &e2t_device, &e12u_device,
            &gphiu_device,
            &rdt, &cbfr, &visc);

        /* Set OpenCL Kernel Parameters for Momentum-v */
        set_args_momv(clkernel[K_MOM_V],
            &width,
            &va_device, &un_device, &vn_device,
            &hu_device, &hv_device, &ht_device,
            &ssha_v_device, &sshn_device,
            &sshn_u_device, &sshn_v_device,
            &tmask_device,
            &e1v_device, &e1t_device,
            &e2u_device, &e2v_device,
            &e2t_device, &e12v_device,
            &gphiv_device,
            &rdt, &cbfr, &visc);

        /* Set OpenCL Kernel Parameters for bc_ssh */
        set_args_bc_ssh(clkernel[K_BC_SSH],
            &width,
            &istep,
            &ssha_device,
            &tmask_device,
            &rdt);

        /* Set OpenCL Kernel Parameters for bc_solid_v */
        set_args_bc_solid_u(clkernel[K_BC_SOLID_U],
            &width,
            &ua_device,
            &tmask_device);

        /* Set OpenCL Kernel Parameters for bc_solid_v */
        set_args_bc_solid_v(clkernel[K_BC_SOLID_V],
            &width,
            &ua_device,
            &tmask_device);

        /* Set OpenCL Kernel Parameters for next_sshu kernel */
        set_args_next_sshu(clkernel[K_NEXT_SSH_U],
            &width, &sshn_u_device, &sshn_device, &tmask_device,
            &e12t_device, &e12u_device);
      
        /* Set OpenCL Kernel Parameters for next_sshv kernel */
        set_args_next_sshv(clkernel[K_NEXT_SSH_V],
            &width, &sshn_v_device, &sshn_device, &tmask_device,
            &e12t_device, &e12v_device);

        /* Create an array to store the event associated with each write
        to the device */
        cl_event *write_events = (cl_event*)malloc(num_buffers*sizeof(cl_event));
        int buf_idx = 0;
        ret = clEnqueueWriteBuffer(command_queue[0], ssha_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_t, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
#ifndef SINGLE_KERNEL
        ret = clEnqueueWriteBuffer(command_queue[0], ssha_u_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_u, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], ssha_v_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_v, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
#endif
/*        ret = clEnqueueWriteBuffer(command_queue[0], sshn_device, 1, 0,
			     (size_t)buff_size, (void *)sshn, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], sshn_u_device, 1, 0,
			     (size_t)buff_size, (void *)sshn_u, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], sshn_v_device, 1, 0,
			     (size_t)buff_size, (void *)sshn_v, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], hu_device, 1, 0,
			     (size_t)buff_size, (void *)hu, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], hv_device, 1, 0,
			     (size_t)buff_size, (void *)hv, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
#ifndef SINGLE_KERNEL
        ret = clEnqueueWriteBuffer(command_queue[0], ht_device, 1, 0,
			     (size_t)buff_size, (void *)ht, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
#endif
        ret = clEnqueueWriteBuffer(command_queue[0], un_device, 1, 0,
			     (size_t)buff_size, (void *)un, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], vn_device, 1, 0,
			     (size_t)buff_size, (void *)vn, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
#ifndef SINGLE_KERNEL
        ret = clEnqueueWriteBuffer(command_queue[0], ua_device, 1, 0,
			     (size_t)buff_size, (void *)ua, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], va_device, 1, 0,
			     (size_t)buff_size, (void *)va, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], e1u_device, 1, 0,
			     (size_t)buff_size, (void *)e1u, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], e1v_device, 1, 0,
			     (size_t)buff_size, (void *)e1v, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], e1t_device, 1, 0,
			     (size_t)buff_size, (void *)e1t, 0,
			     NULL, &(write_events[buf_idx++]));
        ret = clEnqueueWriteBuffer(command_queue[0], e2u_device, 1, 0,
			     (size_t)buff_size, (void *)e2u, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], e2v_device, 1, 0,
			     (size_t)buff_size, (void *)e2v, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], e2t_device, 1, 0,
			     (size_t)buff_size, (void *)e2t, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], e12u_device, 1, 0,
			     (size_t)buff_size, (void *)e12u, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], e12v_device, 1, 0,
			     (size_t)buff_size, (void *)e12v, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
#endif
        ret = clEnqueueWriteBuffer(command_queue[0], e12t_device, 1, 0,
			     (size_t)buff_size, (void *)e12t, 0,
			     NULL, &(write_events[buf_idx++]));
    check_status("clEnqueueWriteBuffer", ret);
#ifndef SINGLE_KERNEL
        ret = clEnqueueWriteBuffer(command_queue[0], gphiu_device, 1, 0,
			     (size_t)buff_size, (void *)gphiu, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], gphiv_device, 1, 0,
			     (size_t)buff_size, (void *)gphiv, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], tmask_device, 1, 0,
			     (size_t)(width*width*sizeof(cl_int)), (void *)tmask, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
#endif
        */
        ret = clWaitForEvents(num_buffers, write_events);
        check_status("clWaitForEvents", ret);


        printf("OpenCL initialization done\n");
        first_time = 0;
    }

    // Continuity kernel (internal domain)
    for(int jj = internal_ystart; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            continuity_code(ji, jj, width, ssha_t, sshn_t, sshn_u, sshn_v, \
                hu, hv, un, vn, rdt, area_t);
        }
    }

    // Momentum_u kernel (internal domain)
    for(int jj = internal_ystart; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            momentum_u_code(ji, jj, width, ua, un, vn, hu, hv, ht, ssha_u, \
                sshn_t, sshn_u, sshn_v, tmask, dx_u, dx_v, dx_t, dy_u, dy_t, \
                area_u, gphiu, rdt, cbfr, visc, omega, d2r, g);
        }
    }

    // Momentum_v kernel (internal domain)
    for(int jj = internal_ystart; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            momentum_v_code(ji, jj, width, va, un, vn, hu, hv, ht, ssha_v, \
                sshn_t, sshn_u, sshn_v, tmask, dx_v, dx_t, dy_u, dy_v, dy_t, \
                area_v, gphiu, rdt, cbfr, visc, omega, d2r, g);
        }
    }

    // Boundary conditions bc_ssh kernel (internal domain)
    for(int jj = internal_ystart; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            bc_ssh_code(ji, jj, width, istep, ssha_t, tmask, rdt);
        }
    }

    // Boundary conditions bc_solid_u kernel (whole domain but top x boundary)
    for(int jj = internal_ystart - 1; jj <= internal_ystop + 1; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop; ji++){
            bc_solid_u_code(ji, jj, width, ua, tmask);
        }
    }

    // Boundary conditions bc_solid_v kernel (whole domain but top y boundary)
    for(int jj = internal_ystart - 1; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop + 1; ji++){
            bc_solid_v_code(ji, jj, width, va, tmask);
        }
    }

    // Boundary conditions bc_flather_u kernel (whole domain but top x boundary)
    for(int jj = internal_ystart - 1; jj <= internal_ystop + 1; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop; ji++){
            bc_flather_u_code(ji, jj, width, ua, hu, sshn_u, tmask, g);
        }
    }

    // Boundary conditions bc_flather_v kernel (whole domain but top y boundary)
    for(int jj = internal_ystart - 1; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop + 1; ji++){
            bc_flather_v_code(ji, jj, width, va, hv, sshn_v, tmask, g);
        }
    }

    // Copy 'next' fields to 'current' fields (whole domain)
    for(int jj = internal_ystart - 1; jj < internal_ystop + 1; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop + 1; ji++){
            int idx = jj * width + ji;
            un[idx] = ua[idx];
            vn[idx] = va[idx];
            sshn_t[idx] = ssha_t[idx];
        }
    }

    // Time update kernel (internal domain u points)
    for(int jj = internal_ystart; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart; ji <= internal_xstop - 1; ji++){
            next_sshu_code(ji, jj, width, sshn_u, sshn_t, tmask, area_t, area_u);
        }
    }

    // Time update kernel (internal domain v points)
    for(int jj = internal_ystart; jj <= internal_ystop - 1; jj++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            next_sshv_code(ji, jj, width, sshn_v, sshn_t, tmask, area_t, area_v);
        }
    }
}
