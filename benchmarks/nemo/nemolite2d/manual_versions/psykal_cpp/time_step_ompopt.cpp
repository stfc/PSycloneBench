#include <iostream>
#include <chrono>
#include <vector>
#include <cstdlib>

// Kernels
#include "../../kernels/c_family/continuity_kern.c"
#include "../../kernels/c_family/momentum_u_kern.c"
#include "../../kernels/c_family/momentum_v_kern.c"
#include "../../kernels/c_family/boundary_conditions_kern.c"
#include "../../kernels/c_family/time_update_kern.c"

// Optimizations:
// - Fuse some kernels (x1.1), fusing more kernels didn't make any difference.
// - Change array copies with swaps is x1.1 faster, but at the moment the
// optimization is commented out because since the fortran->c call passes
// pointers by value the swap is not permanent (produces wrong results).
// - Tests with __restrict__ or assume_aligned didn't improve performance.

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

    // Continuity, momentum and bc_ssh kernels (internal domain)
    #pragma omp parallel for
    for(int jj = internal_ystart; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            continuity_code(ji, jj, width, ssha_t, sshn_t, sshn_u, sshn_v, \
                hu, hv, un, vn, rdt, area_t);
            bc_ssh_code(ji, jj, width, istep, ssha_t, tmask, rdt);
            momentum_u_code(ji, jj, width, ua, un, vn, hu, hv, ht, ssha_u, \
                sshn_t, sshn_u, sshn_v, tmask, dx_u, dx_v, dx_t, dy_u, dy_t, \
                area_u, gphiu, rdt, cbfr, visc, omega, d2r, g);
            momentum_v_code(ji, jj, width, va, un, vn, hu, hv, ht, ssha_v, \
                sshn_t, sshn_u, sshn_v, tmask, dx_v, dx_t, dy_u, dy_v, dy_t, \
                area_v, gphiu, rdt, cbfr, visc, omega, d2r, g);
        }
    }

    // Boundary conditions bc_solid_u kernel (whole domain but top x boundary)
    #pragma omp parallel for
    for(int jj = internal_ystart - 1; jj <= internal_ystop + 1; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop; ji++){
            bc_solid_u_code(ji, jj, width, ua, tmask);
        }
    }

    // Boundary conditions bc_solid_v kernel (whole domain but top y boundary)
    #pragma omp parallel for
    for(int jj = internal_ystart - 1; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop + 1; ji++){
            bc_solid_v_code(ji, jj, width, va, tmask);
        }
    }


    // Boundary conditions bc_flather_u kernel (whole domain but top x boundary)
    #pragma omp parallel for
    for(int jj = internal_ystart - 1; jj <= internal_ystop + 1; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop; ji++){
            bc_flather_u_code(ji, jj, width, ua, hu, sshn_u, tmask, g);
        }
    }

    // Boundary conditions bc_flather_v kernel (whole domain but top y boundary)
    #pragma omp parallel for
    for(int jj = internal_ystart - 1; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop + 1; ji++){
            bc_flather_v_code(ji, jj, width, va, hv, sshn_v, tmask, g);
        }
    }

    // Copy 'next' fields to 'current' fields (whole domain)
    #pragma omp parallel for
    for(int jj = internal_ystart - 1; jj < internal_ystop + 1; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop + 1; ji++){
            int idx = jj * width + ji;
            un[idx] = ua[idx];
            vn[idx] = va[idx];
            sshn_t[idx] = ssha_t[idx];
        }
    }

    // The field copies can be done with a Swap operation as shown below.
    // But at the moment the swapping is not permanent when the
    // function returns to time_step_mod (Fortran PSy layer).
    /*
    double * tmp;
    tmp = un;
    un = ua;
    ua = tmp;

    tmp = vn;
    vn = va;
    va = tmp;

    tmp = sshn_t;
    sshn_t = ssha_t;
    ssha_t = tmp;
    */

    // Time update kernels (internal domain u points)
    #pragma omp parallel for
    for(int jj = internal_ystart; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart; ji <= internal_xstop - 1; ji++){
            next_sshu_code(ji, jj, width, sshn_u, sshn_t, tmask, area_t, area_u);
            next_sshv_code(ji, jj, width, sshn_v, sshn_t, tmask, area_t, area_v);
        }
    }
}
