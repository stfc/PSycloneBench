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
        int total_size,
        double rdt,
        double cbfr,
        double visc,
        double omega,
        double d2r,
        double g
        ){

    // MDRangePolicy uses an open interval (does not include the end
    // point), while the provided 'stop' represent closed ranges.
    // Therefore we need to increase by 1 the 'stop' values (since they
    // are passed by values this only affects this function).
    internal_ystop = internal_ystop + 1;
    internal_xstop = internal_xstop + 1;

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

    // Fields
    Kokkos::View<double**> ssha_t_view("ssha_t", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> sshn_t_view("sshn_t", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> sshn_u_view("sshn_u", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> sshn_v_view("sshn_v", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> hu_view("hu", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> hv_view("hv", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> un_view("un", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> vn_view("vn", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> ua_view("ua", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> ht_view("ht", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> ssha_u_view("ssha_u", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> va_view("va", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> ssha_v_view("ssha_v", internal_xstop+1, internal_ystop+1);

    // Grid
    Kokkos::View<int**> tmask_view("tmask_v", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> area_t_view("area_t", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> area_u_view("area_u", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> area_v_view("area_v", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> dx_u_view("dx_u", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> dx_v_view("dx_v", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> dx_t_view("dx_t", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> dy_u_view("dy_u", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> dy_v_view("dy_v", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> dy_t_view("dy_t", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> gphiu_view("gphiu", internal_xstop+1, internal_ystop+1);
    Kokkos::View<double**> gphiv_view("gphiv", internal_xstop+1, internal_ystop+1);

    // In this implementation the kernels are manually inlined because:
    // - Using the kernels/c_family/* , the GPU version created the LAMBDA and then
    // a call to a gpu function.
    // - It can use 2D Views (which have a (x,y) notation instead of []).

    // Continuity kernel (internal domain)
    Kokkos::parallel_for("continuity",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            double rtmp1, rtmp2, rtmp3, rtmp4;

            rtmp1 = (sshn_u_view(jj, ji)   + hu_view(jj, ji))   * un_view(jj, ji);
            rtmp2 = (sshn_u_view(jj, ji-1) + hu_view(jj, ji-1)) * un_view(jj, ji-1);
            rtmp3 = (sshn_v_view(jj, ji)   + hv_view(jj, ji))   * vn_view(jj, ji);
            rtmp4 = (sshn_v_view(jj-1, ji) + hv_view(jj-1, ji)) * vn_view(jj-1, ji);

            ssha_t_view(jj, ji) = sshn_t_view(jj, ji) + (rtmp2 - rtmp1 + rtmp4 - rtmp3) *
                rdt / area_t_view(jj, ji);
        }
    );
/*
    // Momentum_u kernel (internal domain u points)
    Kokkos::parallel_for("momentum_u",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop - 1}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            momentum_u_code(ji, jj, width, ua_view.data(), un_view.data(), vn_view.data(),
                    hu_view.data(), hv_view.data(), ht_view.data(), ssha_u_view.data(), \
                sshn_t_view.data(), sshn_u_view.data(), sshn_v_view.data(), tmask_view.data(), dx_u_view.data(), dx_v_view.data(), dx_t_view.data(), dy_u_view.data(), dy_t_view.data(), \
                area_u_view.data(), gphiu_view.data(), rdt, cbfr, visc, omega, d2r, g);

        }
    );

    // Momentum_v kernel (internal domain v points)
    Kokkos::parallel_for("momentum_v",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop - 1, internal_xstop}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            momentum_v_code(ji, jj, width, va_view.data(), un_view.data(), vn_view.data(), hu_view.data(), hv_view.data(), ht_view.data(), ssha_v_view.data(), \
                sshn_t_view.data(), sshn_u_view.data(), sshn_v_view.data(), tmask_view.data(), dx_v_view.data(), dx_t_view.data(), dy_u_view.data(), dy_v_view.data(), dy_t_view.data(), \
                area_v_view.data(), gphiv_view.data(), rdt, cbfr, visc, omega, d2r, g);

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
            un(jj, ji) = ua(jj, ji);
            vn(jj, ji) = va(jj, ji);
            sshn_t(jj, ji) = ssha_t(jj, ji);
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
*/
}
