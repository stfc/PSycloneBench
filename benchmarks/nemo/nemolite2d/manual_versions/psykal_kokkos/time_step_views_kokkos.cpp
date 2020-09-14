#include <iostream>
#include <chrono>
#include <vector>
#include <cstdlib>
#include <Kokkos_Core.hpp>

// Uncomment line below to use TIMER
// #define USE_TIMER

#ifdef USE_TIMER
#include "timing.h"
#endif

#define TILE {64,4}

// Create 2D View types for the Fields and Grid arrays
typedef Kokkos::View<double**> double_2dview;
typedef Kokkos::View<int**> int_2dview;

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
        // Device Pointers are passed as void* and the invoke will cast them
        // to the desired type, and by reference as the value has to be stored.
        void * &ssha_t_view_p,
        void * &sshn_t_view_p,
        void * &sshn_u_view_p,
        void * &sshn_v_view_p,
        void * &hu_view_p,
        void * &hv_view_p,
        void * &un_view_p,
        void * &vn_view_p,
        void * &ua_view_p,
        void * &ht_view_p,
        void * &ssha_u_view_p,
        void * &va_view_p,
        void * &ssha_v_view_p,
        void * &tmask_view_p,
        void * &area_t_view_p,
        void * &area_u_view_p,
        void * &area_v_view_p,
        void * &dx_u_view_p,
        void * &dx_v_view_p,
        void * &dx_t_view_p,
        void * &dy_u_view_p,
        void * &dy_v_view_p,
        void * &dy_t_view_p,
        void * &gphiu_view_p,
        void * &gphiv_view_p,
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

#ifdef USE_TIMER
    TimerInit();
    TimerStart("Kokkos SetUp");
#endif

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
    }


    // The execution space is given as a preprocessor define when compiling
    // this file. e.g. `g++ -DEXEC_SPACE=OpenMP time_step_kokkos.cpp -c`
#if defined (EXEC_SPACE)
    using execution_space = Kokkos::EXEC_SPACE;
#else
    using execution_space = Kokkos::DefaultExecutionSpace;
#endif

    if(first_time){
        std::cout << typeid(execution_space).name() << std::endl;
        std::cout << typeid(double_2dview::memory_space).name() << std::endl;

        // Allocate Fields
        ssha_t_view_p = new double_2dview("ssha_t", internal_xstop+1, internal_ystop+1);
        sshn_t_view_p = new double_2dview("sshn_t", internal_xstop+1, internal_ystop+1);
        sshn_u_view_p = new double_2dview("sshn_u", internal_xstop+1, internal_ystop+1);
        sshn_v_view_p = new double_2dview("sshn_v", internal_xstop+1, internal_ystop+1);
        hu_view_p = new double_2dview("hu", internal_xstop+1, internal_ystop+1);
        hv_view_p = new double_2dview("hv", internal_xstop+1, internal_ystop+1);
        un_view_p = new double_2dview("un", internal_xstop+1, internal_ystop+1);
        vn_view_p = new double_2dview("vn", internal_xstop+1, internal_ystop+1);
        ua_view_p = new double_2dview("ua", internal_xstop+1, internal_ystop+1);
        ht_view_p = new double_2dview("ht", internal_xstop+1, internal_ystop+1);
        ssha_u_view_p = new double_2dview("ssha_u", internal_xstop+1, internal_ystop+1);
        va_view_p = new double_2dview("va", internal_xstop+1, internal_ystop+1);
        ssha_v_view_p = new double_2dview("ssha_v", internal_xstop+1, internal_ystop+1);

        // Allocate Grid
        tmask_view_p = new int_2dview("tmask_v", internal_xstop+1, internal_ystop+1);
        area_t_view_p = new double_2dview("area_t", internal_xstop+1, internal_ystop+1);
        area_u_view_p = new double_2dview("area_u", internal_xstop+1, internal_ystop+1);
        area_v_view_p = new double_2dview("area_v", internal_xstop+1, internal_ystop+1);
        dx_u_view_p = new double_2dview("dx_u", internal_xstop+1, internal_ystop+1);
        dx_v_view_p = new double_2dview("dx_v", internal_xstop+1, internal_ystop+1);
        dx_t_view_p = new double_2dview("dx_t", internal_xstop+1, internal_ystop+1);
        dy_u_view_p = new double_2dview("dy_u", internal_xstop+1, internal_ystop+1);
        dy_v_view_p = new double_2dview("dy_v", internal_xstop+1, internal_ystop+1);
        dy_t_view_p = new double_2dview("dy_t", internal_xstop+1, internal_ystop+1);
        gphiu_view_p = new double_2dview("gphiu", internal_xstop+1, internal_ystop+1);
        gphiv_view_p = new double_2dview("gphiv", internal_xstop+1, internal_ystop+1);

    }

    // Create references to avoid using the pointer notation for the views below.
    auto& ssha_t_view =  *static_cast<double_2dview*>(ssha_t_view_p);
    auto& sshn_t_view = *static_cast<double_2dview*>(sshn_t_view_p);
    auto& sshn_u_view = *static_cast<double_2dview*>(sshn_u_view_p);
    auto& sshn_v_view = *static_cast<double_2dview*>(sshn_v_view_p);
    auto& hu_view = *static_cast<double_2dview*>(hu_view_p);
    auto& hv_view = *static_cast<double_2dview*>(hv_view_p);
    auto& un_view = *static_cast<double_2dview*>(un_view_p);
    auto& vn_view = *static_cast<double_2dview*>(vn_view_p);
    auto& ua_view = *static_cast<double_2dview*>(ua_view_p);
    auto& ht_view = *static_cast<double_2dview*>(ht_view_p);
    auto& ssha_u_view = *static_cast<double_2dview*>(ssha_u_view_p);
    auto& va_view = *static_cast<double_2dview*>(va_view_p);
    auto& ssha_v_view = *static_cast<double_2dview*>(ssha_v_view_p);

    auto& tmask_view = *static_cast<int_2dview*>(tmask_view_p);
    auto& area_t_view = *static_cast<double_2dview*>(area_t_view_p);
    auto& area_u_view = *static_cast<double_2dview*>(area_u_view_p);
    auto& area_v_view = *static_cast<double_2dview*>(area_v_view_p);
    auto& dx_u_view = *static_cast<double_2dview*>(dx_u_view_p);
    auto& dx_v_view = *static_cast<double_2dview*>(dx_v_view_p);
    auto& dx_t_view = *static_cast<double_2dview*>(dx_t_view_p);
    auto& dy_u_view = *static_cast<double_2dview*>(dy_u_view_p);
    auto& dy_v_view = *static_cast<double_2dview*>(dy_v_view_p);
    auto& dy_t_view = *static_cast<double_2dview*>(dy_t_view_p);
    auto& gphiu_view = *static_cast<double_2dview*>(gphiu_view_p);
    auto& gphiv_view = *static_cast<double_2dview*>(gphiv_view_p);

    if(first_time){
        // Create Mirrors. These are needed when the execution devices do not
        // share the same memory space as the host, the mirrors synchronise the
        // data in both devices when requested by the deep_copy method. If the
        // execution device has access to the host memory the Mirror overlaps
        // with the View memory location (thus avoiding overheads).
        auto h_ssha_t = Kokkos::create_mirror_view( ssha_t_view );
        auto h_sshn_t = Kokkos::create_mirror_view( sshn_t_view );
        auto h_sshn_u = Kokkos::create_mirror_view( sshn_u_view );
        auto h_sshn_v = Kokkos::create_mirror_view( sshn_v_view );
        auto h_hu = Kokkos::create_mirror_view( hu_view );
        auto h_hv = Kokkos::create_mirror_view( hv_view );
        auto h_un = Kokkos::create_mirror_view( un_view );
        auto h_vn = Kokkos::create_mirror_view( vn_view );
        auto h_ua = Kokkos::create_mirror_view( ua_view );
        auto h_ht = Kokkos::create_mirror_view( ht_view );
        auto h_ssha_u = Kokkos::create_mirror_view( ssha_u_view );
        auto h_va = Kokkos::create_mirror_view( va_view );
        auto h_ssha_v = Kokkos::create_mirror_view( ssha_v_view );

        auto h_tmask = Kokkos::create_mirror_view( tmask_view );
        auto h_area_t = Kokkos::create_mirror_view( area_t_view );
        auto h_area_u = Kokkos::create_mirror_view( area_u_view );
        auto h_area_v = Kokkos::create_mirror_view( area_v_view );
        auto h_dx_u = Kokkos::create_mirror_view( dx_u_view );
        auto h_dx_v = Kokkos::create_mirror_view( dx_v_view );
        auto h_dx_t = Kokkos::create_mirror_view( dx_t_view );
        auto h_dy_u = Kokkos::create_mirror_view( dy_u_view );
        auto h_dy_v = Kokkos::create_mirror_view( dy_v_view );
        auto h_dy_t = Kokkos::create_mirror_view( dy_t_view );
        auto h_gphiu = Kokkos::create_mirror_view( gphiu_view );
        auto h_gphiv = Kokkos::create_mirror_view( gphiv_view );

        // Copy Fortran arrays into the Kokkos View Mirrors
        for(int jj=0; jj < internal_ystop+1; jj++){
            for(int ji=0; ji < internal_xstop+1; ji++){
                int idx = jj*width + ji;
                h_ssha_t(jj, ji) = ssha_t[idx];
                h_sshn_t(jj, ji) = sshn_t[idx];
                h_sshn_u(jj, ji) = sshn_u[idx];
                h_sshn_v(jj, ji) = sshn_v[idx];
                h_hu(jj, ji) = hu[idx];
                h_hv(jj, ji) = hv[idx];
                h_un(jj, ji) = un[idx];
                h_vn(jj, ji) = vn[idx];
                h_ua(jj, ji) = ua[idx];
                h_ht(jj, ji) = ht[idx];
                h_ssha_u(jj, ji) = ssha_u[idx];
                h_va(jj, ji) = va[idx];
                h_ssha_v(jj, ji) = ssha_v[idx];

                h_tmask(jj, ji) = tmask[idx];
                h_area_t(jj, ji) = area_t[idx];
                h_area_u(jj, ji) = area_u[idx];
                h_area_v(jj, ji) = area_v[idx];
                h_dx_u(jj, ji) = dx_u[idx];
                h_dx_v(jj, ji) = dx_v[idx];
                h_dx_t(jj, ji) = dx_t[idx];
                h_dy_u(jj, ji) = dy_u[idx];
                h_dy_v(jj, ji) = dy_v[idx];
                h_dy_t(jj, ji) = dy_t[idx];
                h_gphiu(jj, ji) = gphiu[idx];
                h_gphiv(jj, ji) = gphiv[idx];
            }
        }

        // Update Views with mirror data (only copies if device is not the host)
        Kokkos::deep_copy( ssha_t_view, h_ssha_t );
        Kokkos::deep_copy( sshn_t_view, h_sshn_t );
        Kokkos::deep_copy( sshn_u_view, h_sshn_u );
        Kokkos::deep_copy( sshn_v_view, h_sshn_v );
        Kokkos::deep_copy( hu_view, h_hu );
        Kokkos::deep_copy( hv_view, h_hv );
        Kokkos::deep_copy( un_view, h_un );
        Kokkos::deep_copy( vn_view, h_vn );
        Kokkos::deep_copy( ua_view, h_ua );
        Kokkos::deep_copy( ht_view, h_ht );
        Kokkos::deep_copy( ssha_u_view, h_ssha_u );
        Kokkos::deep_copy( va_view, h_va );
        Kokkos::deep_copy( ssha_v_view, h_ssha_v );

        Kokkos::deep_copy( tmask_view, h_tmask );
        Kokkos::deep_copy( area_t_view, h_area_t );
        Kokkos::deep_copy( area_u_view, h_area_u );
        Kokkos::deep_copy( area_v_view, h_area_v );
        Kokkos::deep_copy( dx_u_view, h_dx_u );
        Kokkos::deep_copy( dx_v_view, h_dx_v );
        Kokkos::deep_copy( dx_t_view, h_dx_t );
        Kokkos::deep_copy( dy_u_view, h_dy_u );
        Kokkos::deep_copy( dy_v_view, h_dy_v );
        Kokkos::deep_copy( dy_t_view, h_dy_t );
        Kokkos::deep_copy( gphiu_view, h_gphiu );
        Kokkos::deep_copy( gphiv_view, h_gphiv );

        first_time = false;
    }

    // In this implementation the kernels are manually inlined because the
    // implementations in ../../kernels/c_family/ use 1D raw pointer syntax
    // instead of 2D Views (which have a (x,y) notation instead of []). These
    // can still be split to a different file using the Kokkos functor pattern
    // if needed for readability.
    
#ifdef USE_TIMER
    TimerStop();
    TimerStart("Continuity Kernel");
#endif

    // Execution policy for a multi-dimensional (2D) iteration space.
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>, execution_space> mdrange_policy;

    // Continuity kernel (internal domain)
    Kokkos::parallel_for("continuity",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop}, TILE),
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


#ifdef USE_TIMER
    TimerStop();
    TimerStart("Momentum Kernels");
#endif

    // Momentum_u kernel (internal domain u points)
    Kokkos::parallel_for("momentum_u",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop - 1}, TILE),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            double u_e, u_w, v_n, v_s;
            double v_nc, v_sc;
            double depe, depw, deps, depn;
            double hpg, adv, cor, vis;
            double dudx_e, dudx_w;
            double dudy_s, dudy_n;
            double uu_e, uu_n, uu_s, uu_w;

            if(tmask_view(jj, ji) + tmask_view(jj, ji + 1) <= 0){
                return; // jump over non-computational domain
            }
            if(tmask_view(jj, ji) <= 0 || tmask_view(jj, ji + 1) <= 0){
                return; // jump over boundary u
            }

            u_e  = 0.5 * (un_view(jj, ji) + un_view(jj, ji + 1)) *
                dy_t_view(jj, ji + 1);   // add length scale.
            depe = ht_view(jj, ji + 1) + sshn_t_view(jj, ji + 1);

            u_w  = 0.5 * (un_view(jj, ji) + un_view(jj, ji - 1)) *
                dy_t_view(jj, ji);   // add length scale
            depw = ht_view(jj, ji) + sshn_t_view(jj, ji);

            v_sc = 0.5 * (vn_view(jj - 1, ji) + vn_view(jj - 1, ji + 1));
            v_s  = 0.5 * v_sc * (dx_v_view(jj - 1, ji) + dx_v_view(jj - 1, ji + 1));
            deps = 0.5 * (hv_view(jj - 1, ji) + sshn_v_view(jj - 1, ji) +
                hv_view(jj - 1, ji + 1) + sshn_v_view(jj - 1, ji + 1));

            v_nc = 0.5 * (vn_view(jj, ji) + vn_view(jj, ji + 1));
            v_n  = 0.5 * v_nc * (dx_v_view(jj, ji) + dx_v_view(jj, ji + 1));
            depn = 0.5 * (hv_view(jj, ji) + sshn_v_view(jj, ji) +
                hv_view(jj, ji + 1) + sshn_v_view(jj, ji + 1));

            // -advection (currently first order upwind)
            uu_w = (0.5 - copysign(0.5, u_w)) * un_view(jj, ji) +
                (0.5 + copysign(0.5, u_w)) * un_view(jj, ji - 1) ;
            uu_e = (0.5 + copysign(0.5, u_e)) * un_view(jj, ji) +
                (0.5 - copysign(0.5, u_e)) * un_view(jj, ji + 1) ;

            if(tmask_view(jj - 1, ji) <=0 || tmask_view(jj - 1, ji + 1) <= 0){
                uu_s = (0.5 - copysign(0.5, v_s)) * un_view(jj, ji);
            }else{
                uu_s = (0.5 - copysign(0.5, v_s)) * un_view(jj, ji) +
                    (0.5 + copysign(0.5, v_s)) * un_view(jj - 1, ji) ;
            }

            if(tmask_view(jj + 1, ji) <=0 || tmask_view(jj + 1 + 1, ji) <= 0){
                uu_n = (0.5 + copysign(0.5, v_n)) * un_view(jj, ji);
            }else{
                uu_n = (0.5 + copysign(0.5, v_n)) * un_view(jj, ji) +
                    (0.5 - copysign(0.5, v_n)) * un_view(jj + 1, ji);
            }

            adv = uu_w * u_w * depw - uu_e * u_e * depe +
                uu_s * v_s * deps - uu_n * v_n * depn;

            // -viscosity

            dudx_e = (un_view(jj, ji + 1) - un_view(jj, ji)) / dx_t_view(jj, ji + 1) *
                (ht_view(jj, ji + 1) + sshn_t_view(jj, ji + 1));
            dudx_w = (un_view(jj, ji) - un_view(jj, ji - 1)) / dx_t_view(jj, ji) *
                (ht_view(jj, ji) + sshn_t_view(jj, ji));
            if(tmask_view(jj - 1, ji) <=0 || tmask_view(jj - 1, ji + 1) <= 0){
                dudy_s = 0.0; // slip boundary
            }else{
                dudy_s = (un_view(jj, ji) - un_view(jj - 1, ji)) / (dy_u_view(jj, ji) +
                    dy_u_view(jj - 1, ji)) * (hu_view(jj, ji) + sshn_u_view(jj, ji) +
                    hu_view(jj - 1, ji) + sshn_u_view(jj - 1, ji));
            }

            if(tmask_view(jj + 1, ji) <= 0 || tmask_view(jj + 1 + 1, ji) <= 0){
                dudy_n = 0.0; // slip boundary
            }else{
                dudy_n = (un_view(jj + 1, ji) - un_view(jj, ji)) / (dy_u_view(jj, ji) +
                        dy_u_view(jj + 1, ji)) * (hu_view(jj, ji) + sshn_u_view(jj, ji) +
                        hu_view(jj + 1, ji) + sshn_u_view(jj + 1, ji));
            }

            vis = (dudx_e - dudx_w ) * dy_u_view(jj, ji)  + 
                (dudy_n - dudy_s ) * dx_u_view(jj, ji) * 0.5;
            vis = visc * vis;   // visc will be an array visc(1:jpijglou) 
                                // for variable viscosity, such as turbulent viscosity

            // -Coriolis' force (can be implemented implicitly)
            cor = 0.5 * (2. * omega * sin(gphiu_view(jj, ji) * d2r) * (v_sc + v_nc)) * 
                area_u_view(jj, ji) * (hu_view(jj, ji) + sshn_u_view(jj, ji));

            // -pressure gradient
            hpg = -g * (hu_view(jj, ji) + sshn_u_view(jj, ji)) * dy_u_view(jj, ji) * 
                (sshn_t_view(jj, ji + 1) - sshn_t_view(jj, ji));

            // -linear bottom friction (implemented implicitly.
            ua_view(jj, ji) = (un_view(jj, ji) * (hu_view(jj, ji) + sshn_u_view(jj, ji)) +
                rdt * (adv + vis + cor + hpg) / area_u_view(jj, ji)) / (hu_view(jj, ji) +
                ssha_u_view(jj, ji)) / (1.0 + cbfr * rdt) ;
        }
    );

    // Momentum_v kernel (internal domain v points)
    Kokkos::parallel_for("momentum_v",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop - 1, internal_xstop}, TILE),
        KOKKOS_LAMBDA (const int jj, const int ji) {

            double u_e, u_w, v_n, v_s;
            double u_ec, u_wc, vv_e, vv_n, vv_s, vv_w;
            double depe, depw, deps, depn;
            double hpg, adv, cor, vis;
            double dvdx_e, dvdx_w, dvdy_n, dvdy_s;
            if(tmask_view(jj, ji) + tmask_view(jj, ji + 1) <= 0){
                return; // jump over non-computatinal domain
            }
            if(tmask_view(jj, ji) <= 0 || tmask_view(jj + 1, ji) <= 0){
                return; // jump over v boundary cells
            }

            // kernel v adv 
            v_n = 0.5 * (vn_view(jj, ji) + vn_view(jj + 1, ji)) *
                dx_t_view(jj + 1, ji); // add length scale.
            depn = ht_view(jj + 1, ji) + sshn_t_view(jj + 1, ji);

            v_s = 0.5 * (vn_view(jj, ji) + vn_view(jj - 1, ji)) *
                dx_t_view(jj, ji); // add length scale
            deps = ht_view(jj, ji) + sshn_t_view(jj, ji);

            u_wc = 0.5 * (un_view(jj, ji - 1) + un_view(jj + 1, ji - 1));
            u_w  = 0.5 * u_wc * (dy_u_view(jj, ji - 1) + dy_u_view(jj + 1, ji - 1));
            depw = 0.50 * (hu_view(jj, ji - 1) + sshn_u_view(jj, ji - 1) + 
                hu_view(jj + 1, ji - 1) + sshn_u_view(jj + 1, ji - 1));

            u_ec = 0.5 * (un_view(jj, ji) + un_view(jj + 1, ji));
            u_e  = 0.5 * u_ec * (dy_u_view(jj, ji) + dy_u_view(jj + 1, ji));
            depe = 0.50 * (hu_view(jj, ji) + sshn_u_view(jj, ji) + 
            hu_view(jj + 1, ji) + sshn_u_view(jj + 1, ji));

            // -advection (currently first order upwind)
            vv_s = (0.5 - copysign(0.5, v_s)) * vn_view(jj, ji) +  
                (0.5 + copysign(0.5, v_s)) * vn_view(jj - 1, ji) ;
            vv_n = (0.5 + copysign(0.5, v_n)) * vn_view(jj, ji) +  
                (0.5 - copysign(0.5, v_n)) * vn_view(jj + 1, ji) ;

            if(tmask_view(jj, ji - 1) <= 0 || tmask_view(jj + 1, ji - 1) <= 0){
                vv_w = (0.5 - copysign(0.5, u_w)) * vn_view(jj, ji);
            }else{
                vv_w = (0.5 - copysign(0.5, u_w)) * vn_view(jj, ji) +  
                    (0.5 + copysign(0.5, u_w)) * vn_view(jj, ji - 1) ;
            }

            if(tmask_view(jj, ji + 1) <= 0 || tmask_view(jj + 1 + 1, ji) <= 0){
                vv_e = (0.5 + copysign(0.5, u_e)) * vn_view(jj, ji);
            }else{
                vv_e = (0.5 + copysign(0.5, u_e)) * vn_view(jj, ji) +
                    (0.5 - copysign(0.5, u_e)) * vn_view(jj, ji + 1);
            }

            adv = vv_w * u_w * depw - vv_e * u_e * depe + 
                vv_s * v_s * deps - vv_n * v_n * depn;

            // -viscosity
            dvdy_n = (vn_view(jj + 1, ji) - vn_view(jj, ji)) / dy_t_view(jj + 1, ji) * 
                (ht_view(jj + 1, ji) + sshn_t_view(jj + 1, ji));
            dvdy_s = (vn_view(jj, ji) - vn_view(jj - 1, ji)) / dy_t_view(jj, ji) * 
                (ht_view(jj, ji) + sshn_t_view(jj, ji));

            if(tmask_view(jj, ji - 1) <= 0 || tmask_view(jj + 1, ji - 1) <= 0){
                dvdx_w = 0.0; // slip boundary
            }else{
                dvdx_w = (vn_view(jj, ji) - vn_view(jj, ji - 1)) / (dx_v_view(jj, ji) +
                    dx_v_view(jj, ji - 1)) * (hv_view(jj, ji) + sshn_v_view(jj, ji) +
                    hv_view(jj, ji - 1) + sshn_v_view(jj, ji - 1));
            }

            if(tmask_view(jj, ji + 1) <= 0 || tmask_view(jj + 1 + 1, ji) <= 0){
                dvdx_e = 0.0; // slip boundary
            }else{
                dvdx_e = (vn_view(jj, ji + 1) - vn_view(jj, ji)) / (dx_v_view(jj, ji) +
                    dx_v_view(jj, ji + 1)) * (hv_view(jj, ji) + sshn_v_view(jj, ji) +
                    hv_view(jj, ji + 1) + sshn_v_view(jj, ji + 1));
            }

            vis = (dvdy_n - dvdy_s ) * dx_v_view(jj, ji)  + 
                (dvdx_e - dvdx_w ) * dy_v_view(jj, ji) * 0.5  ;

            vis = visc * vis;   // visc will be a array visc(1:jpijglou) 
                                // for variable viscosity, such as turbulent viscosity

            // -Coriolis' force (can be implemented implicitly)
            cor = -0.5*(2. * omega * sin(gphiv_view(jj, ji) * d2r) * (u_ec + u_wc)) * 
                area_v_view(jj, ji) * (hv_view(jj, ji) + sshn_v_view(jj, ji));

            // -pressure gradient
            hpg = -g * (hv_view(jj, ji) + sshn_v_view(jj, ji)) * dx_v_view(jj, ji) * 
                (sshn_t_view(jj + 1, ji) - sshn_t_view(jj, ji));

            // -linear bottom friction (implemented implicitly.
            va_view(jj, ji) = (vn_view(jj, ji) * (hv_view(jj, ji) + sshn_v_view(jj, ji)) + 
                rdt * (adv + vis + cor + hpg) / area_v_view(jj, ji) ) / 
                ((hv_view(jj, ji) + ssha_v_view(jj, ji))) / (1.0 + cbfr * rdt) ;
            }
    );
    
#ifdef USE_TIMER
    TimerStop();
    TimerStart("Remaining Kernels");
#endif

    // Boundary conditions bc_ssh kernel (internal domain)
    Kokkos::parallel_for("bc_ssh",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop}, TILE),
        KOKKOS_LAMBDA (const int jj, const int ji) {

            double amp_tide, omega_tide, rtime;

            amp_tide = 0.2;
            omega_tide = 2.0 * 3.14159 / (12.42 * 3600.0);
            rtime = istep * rdt;
  
            if(tmask_view(jj, ji) <= 0) return;
  
            if(tmask_view(jj-1, ji) < 0){
                ssha_t_view(jj, ji) = amp_tide * sin(omega_tide * rtime);
            }else if(tmask_view(jj+1, ji) < 0){
                ssha_t_view(jj, ji) = amp_tide * sin(omega_tide * rtime);
            }else if(tmask_view(jj, ji+1) < 0){
                ssha_t_view(jj, ji) = amp_tide * sin(omega_tide * rtime);
            }else if(tmask_view(jj, ji-1) < 0){
                ssha_t_view(jj, ji) = amp_tide * sin(omega_tide * rtime);
            }
        }
    );
    
    // Boundary conditions bc_solid_u kernel (whole domain but top x boundary)
    Kokkos::parallel_for("bc_solid_u",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop + 1, internal_xstop}, TILE),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            if(tmask_view(jj, ji) * tmask_view(jj, ji+1) == 0){
                ua_view(jj, ji) = 0.0;
            }
        }
    );

    // Boundary conditions bc_solid_v kernel (whole domain but top y boundary)
    Kokkos::parallel_for("bc_solid_v",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop, internal_xstop + 1}, TILE),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            if(tmask_view(jj, ji) * tmask_view(jj+1, ji) == 0){
                va_view(jj, ji) = 0.0;
            }
        }
    );

    // Boundary conditions bc_flather_u kernel (whole domain but top x boundary)
    Kokkos::parallel_for("bc_flather_u",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop + 1, internal_xstop}, TILE),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            //                                   Du                 Dssh
            // Flather open boundary condition [---- = sqrt(g/H) * ------]
            //                                   Dn                 Dn
            // ua and va in du/dn should be the specified tidal forcing

            // Check whether this point lies within the domain
            if(tmask_view(jj, ji) + tmask_view(jj, ji+1) <= -1) return;

            if(tmask_view(jj, ji) < 0){
                ua_view(jj, ji) = ua_view(jj, ji+1) + sqrt(g/hu_view(jj, ji)) *
                    (sshn_u_view(jj, ji) - sshn_u_view(jj, ji+1));
            }else if(tmask_view(jj, ji+1) < 0){
                ua_view(jj, ji) = ua_view(jj, ji-1) + sqrt(g/hu_view(jj, ji)) *
                    (sshn_u_view(jj, ji) - sshn_u_view(jj, ji-1));
            }
        }
    );

    // Boundary conditions bc_flather_v kernel (whole domain but top y boundary)
    Kokkos::parallel_for("bc_flather_v",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop, internal_xstop + 1}, TILE),
        KOKKOS_LAMBDA (const int jj, const int ji) {

            // Check whether this point is inside the simulated domain
            // \todo I could set-up a V-mask using exactly the same code structure
            // as below. Could then apply the BC and multiply by V-mask and thus
            // remove conditionals => get vectorisation.*/
            if(tmask_view(jj, ji) + tmask_view(jj+1, ji) <= -1) return;
    
            if(tmask_view(jj, ji) < 0){
                va_view(jj, ji) = va_view(jj+1, ji) + sqrt(g/hv_view(jj, ji)) *
                    (sshn_v_view(jj, ji) - sshn_v_view(jj+1, ji));
            }else if(tmask_view(jj+1, ji) < 0){
                va_view(jj, ji) = va_view(jj-1, ji) + sqrt(g/hv_view(jj, ji)) *
                    (sshn_v_view(jj, ji) - sshn_v_view(jj-1, ji));
            }
        }
    );

    // Copy 'next' fields to 'current' fields (whole domain)
    Kokkos::parallel_for("copy_fields",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop + 1, internal_xstop + 1}, TILE),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            un_view(jj, ji) = ua_view(jj, ji);
            vn_view(jj, ji) = va_view (jj, ji);
            sshn_t_view(jj, ji) = ssha_t_view(jj, ji);
        }
    );

    // Time update kernel (internal domain u points)
    Kokkos::parallel_for("next_sshu",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop - 1}, TILE),
        KOKKOS_LAMBDA (const int jj, const int ji) {

            double rtmp1;

            if(tmask_view(jj, ji) + tmask_view(jj, ji + 1) <= 0){
                return; // jump over non-computational domain
            }

            if(tmask_view(jj, ji) * tmask_view(jj, ji + 1) > 0){
                rtmp1 = area_t_view(jj, ji) * sshn_t_view(jj, ji) +
                    area_t_view(jj, ji + 1) * sshn_t_view(jj, ji + 1);
                sshn_u_view(jj, ji) = 0.5 * rtmp1 / area_u_view(jj, ji) ;
            }else if(tmask_view(jj, ji) <= 0){
                sshn_u_view(jj, ji) = sshn_t_view(jj, ji + 1);
            }else if(tmask_view(jj, ji + 1) <= 0){
                sshn_u_view(jj, ji) = sshn_t_view(jj, ji);
            }
        }
    );

    // Time update kernel (internal domain v points)
    Kokkos::parallel_for("next_sshv",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop - 1, internal_xstop}, TILE),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            double rtmp1;
            if((tmask_view(jj, ji) + tmask_view(jj + 1, ji)) <= 0){
                return; //jump over non-computational domain
            }
            if((tmask_view(jj, ji) * tmask_view(jj + 1, ji)) > 0){
                rtmp1 = area_t_view(jj, ji) * sshn_t_view(jj, ji) +
                    area_t_view(jj + 1, ji) * sshn_t_view(jj + 1, ji);
                sshn_v_view(jj, ji) = 0.5 * rtmp1 / area_v_view(jj, ji) ;
            }else if(tmask_view(jj, ji) <= 0){
                sshn_v_view(jj, ji) = sshn_t_view(jj + 1, ji);
            }else if(tmask_view(jj + 1, ji) <= 0){
                sshn_v_view(jj, ji) = sshn_t_view(jj, ji);
            }
        }
    );

#ifdef USE_TIMER
    TimerStop();
    TimerReport();
#endif

}

extern "C" void kokkos_read_from_device(double_2dview from, double* to,
                                        int nx, int ny, int width){

    // First, we need to copy data from the device into a host mirror ( which
    // has the same data-layout as the device copy), it requires to allocate
    // space for the mirror with `create_mirror_view`
    auto mirror = Kokkos::create_mirror_view( from );
    Kokkos::deep_copy( mirror, from );

    // Then, we copy the data from the mirror to the original location.
    // Since the mirror data layout is decided by kokkos, we make explicit
    // copies of each element to its location.
    for(int jj=0; jj < ny; jj++){
        for(int ji=0; ji < nx; ji++){
            int idx = jj*width + ji;
            to[idx] = mirror(jj, ji);
        }
    }
}
