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
cl_mem area_t_device = NULL;
cl_mem area_u_device = NULL;
cl_mem area_v_device = NULL;
cl_mem dx_u_device = NULL;
cl_mem dx_v_device = NULL;
cl_mem dx_t_device = NULL;
cl_mem dy_u_device = NULL;
cl_mem dy_v_device = NULL;
cl_mem dy_t_device = NULL;
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
        double * gphiv,
        // Scalars
        int istep,
        int internal_xstart,
        int internal_xstop,
        int internal_ystart,
        int internal_ystop,
        int width,
        int total_size,
        double rdt,
        double cbfr,
        double visc,
        double omega,
        double d2r,
        double g
        ){

    int ret;
    int buff_size = total_size * sizeof(cl_double);

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
  
        /* Create Device Memory Buffers*/
        int num_buffers = 0;
        fprintf(stdout, "Creating buffers of size %d ...\n", total_size);
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
        ua_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        va_device = clCreateBuffer(context, CL_MEM_READ_WRITE, buff_size,
			                       NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);

        /* Mesh scale factors */
        area_t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        area_u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        area_v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        dx_u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        dx_v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        dx_t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
			                        NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        dy_u_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                                     NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        dy_v_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                                     NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        dy_t_device = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size,
                                     NULL, &ret);
        num_buffers++;
        check_status("clCreateBuffer", ret);
        tmask_device = clCreateBuffer(context, CL_MEM_READ_ONLY,
				        (size_t)(total_size*sizeof(cl_int)), NULL, &ret);
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
		    &rdt, &area_t_device); 
        /* Set OpenCL Kernel Parameters for Momentum-u */
        set_args_momu(clkernel[K_MOM_U],
            &width,
            &ua_device, &un_device, &vn_device,
            &hu_device, &hv_device, &ht_device,
            &ssha_u_device, &sshn_device,
            &sshn_u_device, &sshn_v_device,
            &tmask_device,
            &dx_u_device, &dx_v_device,
            &dx_t_device, &dy_u_device,
            &dy_t_device, &area_u_device,
            &gphiu_device,
            &rdt, &cbfr, &visc,
            &omega, &d2r, &g);

        /* Set OpenCL Kernel Parameters for Momentum-v */
        set_args_momv(clkernel[K_MOM_V],
            &width,
            &va_device, &un_device, &vn_device,
            &hu_device, &hv_device, &ht_device,
            &ssha_v_device, &sshn_device,
            &sshn_u_device, &sshn_v_device,
            &tmask_device,
            &dx_v_device, &dx_t_device,
            &dy_u_device, &dy_v_device,
            &dy_t_device, &area_v_device,
            &gphiv_device,
            &rdt, &cbfr, &visc,
            &omega, &d2r, &g);

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

        /* Set OpenCL Kernel Parameters for bc_flather_u */
        set_args_bc_flather_u(clkernel[K_BC_FLATHER_U],
            &width, &ua_device, &hu_device,
            &sshn_u_device, &tmask_device, &g);

        /* Set OpenCL Kernel Parameters for bc_flather_v */
        set_args_bc_flather_v(clkernel[K_BC_FLATHER_V],
            &width, &va_device, &hv_device,
            &sshn_v_device, &tmask_device, &g);

        /* Set OpenCL Kernel Parameters for next_sshu kernel */
        set_args_next_sshu(clkernel[K_NEXT_SSH_U],
            &width, &sshn_u_device, &sshn_device, &tmask_device,
            &area_t_device, &area_u_device);
      
        /* Set OpenCL Kernel Parameters for next_sshv kernel */
        set_args_next_sshv(clkernel[K_NEXT_SSH_V],
            &width, &sshn_v_device, &sshn_device, &tmask_device,
            &area_t_device, &area_v_device);

        /* Create an array to store the event associated with each write
        to the device */
        cl_event *write_events = (cl_event*)malloc(num_buffers*sizeof(cl_event));
        int buf_idx = 0;
        // Fields
        ret = clEnqueueWriteBuffer(command_queue[0], ssha_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_t, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], ssha_u_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_u, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], ssha_v_device, 1, 0,
			     (size_t)buff_size, (void *)ssha_v, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], sshn_device, 1, 0,
			     (size_t)buff_size, (void *)sshn_t, 0,
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
        ret = clEnqueueWriteBuffer(command_queue[0], ht_device, 1, 0,
			     (size_t)buff_size, (void *)ht, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], un_device, 1, 0,
			     (size_t)buff_size, (void *)un, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], vn_device, 1, 0,
			     (size_t)buff_size, (void *)vn, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], ua_device, 1, 0,
			     (size_t)buff_size, (void *)ua, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], va_device, 1, 0,
			     (size_t)buff_size, (void *)va, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);

        // Grid
        ret = clEnqueueWriteBuffer(command_queue[0], area_t_device, 1, 0,
			     (size_t)buff_size, (void *)area_t, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], area_u_device, 1, 0,
			     (size_t)buff_size, (void *)area_u, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], area_v_device, 1, 0,
			     (size_t)buff_size, (void *)area_v, 0,
			     NULL, &(write_events[buf_idx++]));
        ret = clEnqueueWriteBuffer(command_queue[0], dx_u_device, 1, 0,
			     (size_t)buff_size, (void *)dx_u, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], dx_v_device, 1, 0,
			     (size_t)buff_size, (void *)dx_v, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], dx_t_device, 1, 0,
			     (size_t)buff_size, (void *)dx_t, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], dy_u_device, 1, 0,
			     (size_t)buff_size, (void *)dy_u, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], dy_v_device, 1, 0,
			     (size_t)buff_size, (void *)dy_v, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], dy_t_device, 1, 0,
			     (size_t)buff_size, (void *)dy_t, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], gphiu_device, 1, 0,
			     (size_t)buff_size, (void *)gphiu, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], gphiv_device, 1, 0,
			     (size_t)buff_size, (void *)gphiv, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clEnqueueWriteBuffer(command_queue[0], tmask_device, 1, 0,
			     (size_t)(total_size*sizeof(cl_int)), (void *)tmask, 0,
			     NULL, &(write_events[buf_idx++]));
        check_status("clEnqueueWriteBuffer", ret);
        ret = clWaitForEvents(buf_idx, write_events);
        check_status("clWaitForEvents", ret);


        printf("OpenCL initialization done\n");
        printf("OpenCL ystart=%d\n",internal_ystart);
        printf("OpenCL ystop=%d\n",internal_ystop);
        printf("OpenCL xstart=%d\n",internal_xstart);
        printf("OpenCL xstop=%d\n",internal_xstop);
        first_time = 0;
    }

    // To get good performance the first dimension of the local size need to be
    // {32, 1} or {64,1}. The current implementation expects the second
    // dimension to always be 1.
    size_t local_size[2] = {64, 1};

    // The work_group_size[0] values below should always be evently divisible by
    // the local_size[0]. For this reason, sometimes it is necessary to increase
    // the xstop value until it satisfies the condition. A new computed_xstop
    // value is computed here which adds between 1 and total_size[0] to xstop.
    // Note that the kernels will need a condition e.g:
    //     if (ji > internal_xstop) return;
    // to mask out the additional iterations.
    int extra_xstop =  local_size[0] - (internal_xstop % local_size[0]);
    int computed_xstop = internal_xstop + extra_xstop;

    // Continuity kernel (internal domain)
    size_t continuity_offset[2] = {(size_t)internal_xstart, (size_t)internal_ystart};
    size_t continuity_size[2] = {(size_t)computed_xstop, (size_t)internal_ystop}; 
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_CONTINUITY], 2,
        continuity_offset, continuity_size, local_size, 0, NULL,
        &(clkernevt[K_CONTINUITY]));
    check_status("clEnqueueNDRangeKernel(Continuity)", ret);

    // Momentum_u kernel (internal domain)
    size_t momu_offset[2] = {(size_t)internal_xstart, (size_t)internal_ystart};
    //size_t momu_size[2] = {(size_t)internal_xstop - 1, (size_t)internal_ystop};
    size_t momu_size[2] = {(size_t)computed_xstop, (size_t)internal_ystop};
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_MOM_U], 2,
            momu_offset, momu_size, local_size, 0, NULL,
			&(clkernevt[K_MOM_U]));
    check_status("clEnqueueNDRangeKernel(Mom-u)", ret);

    // Momentum_v kernel (internal domain)
    size_t momv_offset[2] = {(size_t)internal_xstart, (size_t)internal_ystart};
    size_t momv_size[2] = {(size_t)computed_xstop, (size_t)internal_ystop - 1};
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_MOM_V], 2,
            momv_offset, momv_size, local_size, 0, NULL,
			&(clkernevt[K_MOM_V]));
    check_status("clEnqueueNDRangeKernel(Mom-v)", ret);

    /* Set bc_ssh_code again because the istep argument is different */
    set_args_bc_ssh(clkernel[K_BC_SSH],
        &width,
        &istep,
        &ssha_device,
        &tmask_device,
        &rdt);

    // Boundary conditions bc_ssh kernel (internal domain)
    size_t bcssh_offset[2] = {(size_t)internal_xstart, (size_t)internal_ystart};
    size_t bcssh_size[2] = {(size_t)computed_xstop, (size_t)internal_ystop}; 
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_BC_SSH], 2,
            bcssh_offset, bcssh_size, local_size, 0, NULL,
			&(clkernevt[K_BC_SSH]));
    check_status("clEnqueueNDRangeKernel(bc-ssh)", ret);

    // Boundary conditions bc_solid_u kernel (whole domain but top x boundary)
    size_t solidu_offset[2] = {(size_t)internal_xstart - 1, (size_t)internal_ystart - 1};
    size_t solidu_size[2] = {(size_t)computed_xstop, (size_t)internal_ystop + 1}; 
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_BC_SOLID_U], 2,
            solidu_offset, solidu_size, local_size, 0, NULL,
			&(clkernevt[K_BC_SOLID_U]));
    check_status("clEnqueueNDRangeKernel(bc-solid-u)", ret);

    // Boundary conditions bc_solid_v kernel (whole domain but top y boundary)
    size_t solidv_offset[2] = {(size_t)internal_xstart - 1, (size_t)internal_ystart - 1};
    //size_t solidv_size[2] = {(size_t)internal_xstop + 1, (size_t)internal_ystop}; 
    size_t solidv_size[2] = {(size_t)computed_xstop, (size_t)internal_ystop}; 
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_BC_SOLID_V], 2,
            solidv_offset, solidv_size, local_size, 0, NULL,
			&(clkernevt[K_BC_SOLID_V]));
    check_status("clEnqueueNDRangeKernel(bc-solid-v)", ret);

    // Boundary conditions bc_flather_u kernel (whole domain but top x boundary)
    size_t flatheru_offset[2] = {(size_t)internal_xstart - 1, (size_t)internal_ystart - 1};
    size_t flatheru_size[2] = {(size_t)computed_xstop, (size_t)internal_ystop + 1}; 
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_BC_FLATHER_U], 2,
            flatheru_offset, flatheru_size, local_size, 0, NULL,
			&(clkernevt[K_BC_FLATHER_U]));
    check_status("clEnqueueNDRangeKernel(bc-flather-u)", ret);

    // Boundary conditions bc_flather_v kernel (whole domain but top y boundary)
    size_t flatherv_offset[2] = {(size_t)internal_xstart - 1, (size_t)internal_ystart - 1};
    //size_t flatherv_size[2] = {(size_t)internal_xstop + 1, (size_t)internal_ystop}; 
    size_t flatherv_size[2] = {(size_t)computed_xstop, (size_t)internal_ystop}; 
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_BC_FLATHER_V], 2,
            flatherv_offset, flatherv_size, local_size, 0, NULL,
			&(clkernevt[K_BC_FLATHER_V]));
    check_status("clEnqueueNDRangeKernel(bc-flather-v)", ret);

    // Copy 'next' fields to 'current' fields (whole domain)
    ret = clEnqueueCopyBuffer(command_queue[0], ua_device, un_device, 0, 0,
			      buff_size,0, NULL, NULL);
    check_status("clEnqueueCopyBuffer", ret);
    ret = clEnqueueCopyBuffer(command_queue[0], va_device, vn_device, 0, 0,
			      buff_size,0, NULL, NULL);
    check_status("clEnqueueCopyBuffer", ret);
    ret = clEnqueueCopyBuffer(command_queue[0], ssha_device, sshn_device, 0, 0,
			      buff_size,0, NULL, NULL);
    check_status("clEnqueueCopyBuffer", ret);

    // Time update kernel (internal domain u points)
    size_t nextu_offset[2] = {(size_t)internal_xstart, (size_t)internal_ystart};
    //size_t nextu_size[2] = {(size_t)internal_xstop - 1, (size_t)internal_ystop}; 
    size_t nextu_size[2] = {(size_t)computed_xstop, (size_t)internal_ystop}; 
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_NEXT_SSH_U], 2,
            nextu_offset, nextu_size, local_size, 0, NULL,
			&(clkernevt[K_NEXT_SSH_U]));
    check_status("clEnqueueNDRangeKernel(next-sshu)", ret);

    // Time update kernel (internal domain v points)
    size_t nextv_offset[2] = {(size_t)internal_xstart, (size_t)internal_ystart};
    size_t nextv_size[2] = {(size_t)computed_xstop, (size_t)internal_ystop - 1};
    ret = clEnqueueNDRangeKernel(command_queue[0], clkernel[K_NEXT_SSH_V], 2,
            nextv_offset, nextv_size, local_size, 0, NULL,
			&(clkernevt[K_NEXT_SSH_V]));
    check_status("clEnqueueNDRangeKernel(next-sshv)", ret);

    // Wait for the completion of all kernels  -- wait below assumes they are in-order
    ret = clWaitForEvents(1, &(clkernevt[K_NEXT_SSH_V]));
    check_status("clWaitForEvents", ret);

    // Read the data back
    cl_event read_events[3];
    int nread = 0;
    ret = clEnqueueReadBuffer(command_queue[0], ssha_device, CL_TRUE, 0,
			    (size_t)buff_size, (void *)ssha_t, 0, NULL, &(read_events[0]));
    nread++;
    check_status("clEnqueueReadBuffer", ret);
    ret = clEnqueueReadBuffer(command_queue[0], ua_device, CL_TRUE, 0,
			    (size_t)buff_size, (void *)ua, 0, NULL, &(read_events[1]));
    nread++;
    check_status("clEnqueueReadBuffer", ret);
    ret = clEnqueueReadBuffer(command_queue[0], va_device, CL_TRUE, 0,
			    (size_t)buff_size, (void *)va, 0, NULL, &(read_events[2]));
    nread++;
    check_status("clEnqueueReadBuffer", ret);

    clWaitForEvents(nread, read_events);
    check_status("clWaitForEvents", ret);

}
