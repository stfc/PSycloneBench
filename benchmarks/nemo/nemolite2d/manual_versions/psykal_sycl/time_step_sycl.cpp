#include <iostream>
#include <chrono>
#include <vector>
#include <cstdlib>
#include <CL/sycl.hpp>

using namespace cl::sycl;

// Uncomment line below to use TIMER
// #define USE_TIMER

#ifdef USE_TIMER
#include "timing.h"
#endif

queue * workqueue;
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
    TimerStart("Sycl SetUp");
#endif

    if(first_time){
        std::cout << "List of SYCL devices:" << std::endl;
        auto devices = device::get_devices();
        for(auto dev : devices){
            auto platform_name = dev.get_platform().get_info<info::platform::name>();
            auto dev_name = dev.get_info<info::device::name>();
            std::cout << " - " << platform_name << " : " << dev_name << std::endl;
        }
        workqueue = new queue(devices[3]);
        first_time = false;
    }

    auto& myqueue = *workqueue;

    // MDRangePolicy uses an open interval (does not include the end
    // point), while the provided 'stop' represent closed ranges.
    // Therefore we need to increase by 1 the 'stop' values (since they
    // are passed by values this only affects this function).
    internal_ystop = internal_ystop + 1;
    internal_xstop = internal_xstop + 1;


    int height = width;
    buffer<double, 2> ssha_t_buffer(ssha_t, range<2>(width,height));
    buffer<double, 2> sshn_t_buffer(sshn_t, range<2>(width, height));
    buffer<double, 2> sshn_u_buffer(sshn_u, range<2>(width, height));
    buffer<double, 2> sshn_v_buffer(sshn_v, range<2>(width, height));
    buffer<double, 2> hu_buffer(hu, range<2>(width, height));
    buffer<double, 2> hv_buffer(hv, range<2>(width, height));
    buffer<double, 2> un_buffer(un, range<2>(width, height));
    buffer<double, 2> vn_buffer(vn, range<2>(width, height));
    buffer<double, 2> area_t_buffer(area_t, range<2>(width, height));
    std::cout << "Cont: " << width << " " << height << std::endl; 
    std::cout << "Cont: " << internal_ystop << " " << internal_xstop << std::endl; 

    std::cout << "Cont: " << ssha_t[0] << " " << ssha_t[1] << std::endl; 
    myqueue.submit([&](handler &cgh){
        auto ssha_t_accessor = ssha_t_buffer.get_access<access::mode::write>(cgh);
        auto sshn_u_accessor = sshn_u_buffer.get_access<access::mode::read>(cgh);
        auto sshn_v_accessor = sshn_v_buffer.get_access<access::mode::read>(cgh);
        auto sshn_t_accessor = sshn_t_buffer.get_access<access::mode::read>(cgh);
        auto hu_accessor = hu_buffer.get_access<access::mode::read>(cgh);
        auto hv_accessor = hv_buffer.get_access<access::mode::read>(cgh);
        auto un_accessor = un_buffer.get_access<access::mode::read>(cgh);
        auto vn_accessor = vn_buffer.get_access<access::mode::read>(cgh);
        auto area_t_accessor = area_t_buffer.get_access<access::mode::read>(cgh);

        cgh.parallel_for(range<2>(internal_ystop, internal_xstop), [=](id<2> idx){
            double rtmp1, rtmp2, rtmp3, rtmp4;
            auto ji = idx[0];
            auto jj = idx[1];
            if (ji < internal_xstart) return;
            if (jj < internal_ystart) return;

            rtmp1 = (sshn_u_accessor[idx]        + hu_accessor[idx])        * un_accessor[idx];
            rtmp2 = (sshn_u_accessor[{jj, ji-1}] + hu_accessor[{jj, ji-1}]) * un_accessor[{jj, ji-1}];
            rtmp3 = (sshn_v_accessor[{jj, ji}]   + hv_accessor[{jj, ji}])   * vn_accessor[{jj, ji}];
            rtmp4 = (sshn_v_accessor[{jj-1, ji}] + hv_accessor[{jj-1, ji}]) * vn_accessor[{jj-1, ji}];

            ssha_t_accessor[idx] = sshn_t_accessor[idx] + (rtmp2 - rtmp1 + rtmp4 - rtmp3) *
                rdt / area_t_accessor[idx];
        });
            
    });

    myqueue.wait();

    //std::cout << "Cont: " << ssha_t[500] << " " << ssha_t[600] << std::endl; 
    /*

    // Execution policy for a multi-dimensional (2D) iteration space.
    typedef Kokkos::MDRangePolicy<Kokkos::Rank<2>, execution_space> mdrange_policy;

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


#ifdef USE_TIMER
    TimerStop();
    TimerStart("Momentum Kernels");
#endif

    // Momentum_u kernel (internal domain u points)
    Kokkos::parallel_for("momentum_u",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop - 1}),
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
                       {internal_ystop - 1, internal_xstop}),
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
                       {internal_ystop, internal_xstop}),
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
                       {internal_ystop + 1, internal_xstop}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            if(tmask_view(jj, ji) * tmask_view(jj, ji+1) == 0){
                ua_view(jj, ji) = 0.0;
            }
        }
    );

    // Boundary conditions bc_solid_v kernel (whole domain but top y boundary)
    Kokkos::parallel_for("bc_solid_v",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop, internal_xstop + 1}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            if(tmask_view(jj, ji) * tmask_view(jj+1, ji) == 0){
                va_view(jj, ji) = 0.0;
            }
        }
    );

    // Boundary conditions bc_flather_u kernel (whole domain but top x boundary)
    Kokkos::parallel_for("bc_flather_u",
        mdrange_policy({internal_ystart - 1, internal_xstart - 1},
                       {internal_ystop + 1, internal_xstop}),
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
                       {internal_ystop, internal_xstop + 1}),
        KOKKOS_LAMBDA (const int jj, const int ji) {

            // Check whether this point is inside the simulated domain
            // \todo I could set-up a V-mask using exactly the same code structure
            // as below. Could then apply the BC and multiply by V-mask and thus
            // remove conditionals => get vectorisation.
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
                       {internal_ystop + 1, internal_xstop + 1}),
        KOKKOS_LAMBDA (const int jj, const int ji) {
            un_view(jj, ji) = ua_view(jj, ji);
            vn_view(jj, ji) = va_view (jj, ji);
            sshn_t_view(jj, ji) = ssha_t_view(jj, ji);
        }
    );

    // Time update kernel (internal domain u points)
    Kokkos::parallel_for("next_sshu",
        mdrange_policy({internal_ystart, internal_xstart},
                       {internal_ystop, internal_xstop - 1}),
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
                       {internal_ystop - 1, internal_xstop}),
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
*/
#ifdef USE_TIMER
    TimerStop();
    TimerReport();
#endif

}

extern "C" void kokkos_read_from_device(void* from, double* to,
                                        int nx, int ny, int width){

}
