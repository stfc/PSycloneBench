#include <iostream>
#include "opencl_utils.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

// Kernels
#include "../../kernels/c_family/continuity_kern.c"
#include "../../kernels/c_family/momentum_u_kern.c"
#include "../../kernels/c_family/momentum_v_kern.c"
#include "../../kernels/c_family/boundary_conditions_kern.c"
#include "../../kernels/c_family/time_update_kern.c"

/** Maximum number of OpenCL devices we will query */
#define MAX_DEVICES 4

/* Number of OpenCL command queues we will use */
#define NUM_QUEUES 2

bool first_time = true;

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



extern "C" void c_invoke_time_step(
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

        init_device(&device, version_str, &context);

        printf("OpenCL initialization done\n");
        first_time = false;
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
