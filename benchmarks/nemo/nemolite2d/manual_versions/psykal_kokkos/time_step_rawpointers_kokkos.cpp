#include <iostream>
#include <chrono>
#include <vector>
#include <cstdlib>
#include <Kokkos_Core.hpp>

// Kernels
#include "../../kernels/c_family/continuity_kern.c"
#include "../../kernels/c_family/momentum_u_kern.c"
#include "../../kernels/c_family/momentum_v_kern.c"
#include "../../kernels/c_family/boundary_conditions_kern.c"
#include "../../kernels/c_family/time_update_kern.c"

bool first_time = true;

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
        double * gphiv,
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

    // Kokkos needs to be initialized. Since NemoLite2D just has a single
    // invoke, it is simple to do it here the first time the invoke is
    // executed. Note that this can not be done for `Kokkos::finalize();`
    // which is ignored in this implementation. 
    if(first_time){
        Kokkos::initialize();
        first_time = false;
    }


    // The execution space is given as a preprocessor define when compiling
    // this file. e.g. `g++ -DEXEC_SPACE=OpenMP time_step_kokkos.cpp -c`
#if defined (EXEC_SPACE)
    using execution_space = Kokkos::EXEC_SPACE;
#else
    using execution_space = Kokkos::DefaultExecutionSpace;
#endif

    // Execution policy for a multi-dimensional (2D) iteration space.
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>, execution_space> mdrange_policy;

    // MDRangePolicy uses an open interval (does not include the end
    // point), while the provided 'stop' represent closed ranges.
    // Therefore we need to increase by 1 the 'stop' values (since they
    // are passed by values this only affects this function).
    internal_ystop = internal_ystop + 1;
    internal_xstop = internal_xstop + 1;

    // Continuity kernel (internal domain)
    Kokkos::parallel_for("continuity",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            continuity_code(ji, jj, width, ssha_t, sshn_t, sshn_u, sshn_v, \
                hu, hv, un, vn, rdt, area_t);
        }
    );

    // Momentum_u kernel (internal domain u points)
    Kokkos::parallel_for("momentum_u",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop - 1}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            momentum_u_code(ji, jj, width, ua, un, vn, hu, hv, ht, ssha_u, \
                sshn_t, sshn_u, sshn_v, tmask, dx_u, dx_v, dx_t, dy_u, dy_t, \
                area_u, gphiu, rdt, cbfr, visc, omega, d2r, g);

        }
    );

    // Momentum_v kernel (internal domain v points)
    Kokkos::parallel_for("momentum_v",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop - 1, internal_xstop}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            momentum_v_code(ji, jj, width, va, un, vn, hu, hv, ht, ssha_v, \
                sshn_t, sshn_u, sshn_v, tmask, dx_v, dx_t, dy_u, dy_v, dy_t, \
                area_v, gphiv, rdt, cbfr, visc, omega, d2r, g);

        }
    );
    
    // Boundary conditions bc_ssh kernel (internal domain)
    Kokkos::parallel_for("bc_ssh",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            bc_ssh_code(ji, jj, width, istep, ssha_t, tmask, rdt);
        }
    );
    
    // Boundary conditions bc_solid_u kernel (whole domain but top x boundary)
    Kokkos::parallel_for("bc_ssh",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop + 1, internal_xstop}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            bc_solid_u_code(ji, jj, width, ua, tmask);
        }
    );

    // Boundary conditions bc_solid_v kernel (whole domain but top y boundary)
    Kokkos::parallel_for("bc_solid_v",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop, internal_xstop + 1}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            bc_solid_u_code(ji, jj, width, ua, tmask);
        }
    );

    // Boundary conditions bc_flather_u kernel (whole domain but top x boundary)
    Kokkos::parallel_for("bc_solid_v",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop + 1, internal_xstop}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            bc_flather_u_code(ji, jj, width, ua, hu, sshn_u, tmask, g);
        }
    );

    // Boundary conditions bc_flather_v kernel (whole domain but top y boundary)
    Kokkos::parallel_for("bc_solid_v",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop, internal_xstop + 1}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            bc_flather_v_code(ji, jj, width, va, hv, sshn_v, tmask, g);
        }
    );

    // Copy 'next' fields to 'current' fields (whole domain)
    Kokkos::parallel_for("copy_fields",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop + 1, internal_xstop + 1}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            int idx = jj * width + ji;
            un[idx] = ua[idx];
            vn[idx] = va[idx];
            sshn_t[idx] = ssha_t[idx];
        }
    );

    // Time update kernel (internal domain u points)
    Kokkos::parallel_for("next_sshu",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop - 1}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            next_sshu_code(ji, jj, width, sshn_u, sshn_t, tmask, area_t, area_u);
        }
    );

    // Time update kernel (internal domain v points)
    Kokkos::parallel_for("next_sshv",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop - 1, internal_xstop}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            next_sshv_code(ji, jj, width, sshn_v, sshn_t, tmask, area_t, area_v);
        }
    );

}
