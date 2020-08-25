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

#include "../../kernels/c_family/momentum_u_kern.c"
#include "../../kernels/c_family/momentum_v_kern.c"

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
        int index = 0;
        for(auto dev : devices){
            auto platform_name = dev.get_platform().get_info<info::platform::name>();
            auto dev_name = dev.get_info<info::device::name>();
            std::cout <<  " " << index++ << " ) " << platform_name << " : " <<
                dev_name << std::endl;
        }

        if(const char* env_p = std::getenv("SYCL_PLATFORM")){
            std::exit(-1);
        } else {
            // Use default device
            workqueue = new queue(devices[3]);
        }
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

    buffer<double, 2> ssha_t_buffer(ssha_t, range<2>(width, height));
    buffer<double, 2> sshn_t_buffer(sshn_t, range<2>(width, height));
    buffer<double, 2> sshn_u_buffer(sshn_u, range<2>(width, height));
    buffer<double, 2> sshn_v_buffer(sshn_v, range<2>(width, height));
    buffer<double, 2> hu_buffer(hu, range<2>(width, height));
    buffer<double, 2> hv_buffer(hv, range<2>(width, height));
    buffer<double, 2> un_buffer(un, range<2>(width, height));
    buffer<double, 2> vn_buffer(vn, range<2>(width, height));
    buffer<double, 2> ua_buffer(ua, range<2>(width, height));
    buffer<double, 2> ht_buffer(ht, range<2>(width, height));
    buffer<double, 2> ssha_u_buffer(ssha_u, range<2>(width, height));
    buffer<double, 2> va_buffer(va, range<2>(width, height));
    buffer<double, 2> ssha_v_buffer(ssha_v, range<2>(width, height));

    buffer<int, 2> tmask_buffer(tmask, range<2>(width, height));
    buffer<double, 2> area_t_buffer(area_t, range<2>(width, height));
    buffer<double, 2> area_u_buffer(area_u, range<2>(width, height));
    buffer<double, 2> area_v_buffer(area_v, range<2>(width, height));
    buffer<double, 2> dx_u_buffer(dx_u, range<2>(width, height));
    buffer<double, 2> dx_v_buffer(dx_v, range<2>(width, height));
    buffer<double, 2> dx_t_buffer(dx_t, range<2>(width, height));
    buffer<double, 2> dy_u_buffer(dy_u, range<2>(width, height));
    buffer<double, 2> dy_v_buffer(dy_v, range<2>(width, height));
    buffer<double, 2> dy_t_buffer(dy_t, range<2>(width, height));
    buffer<double, 2> gphiu_buffer(gphiu, range<2>(width, height));
    buffer<double, 2> gphiv_buffer(gphiv, range<2>(width, height));

    myqueue.submit([&](handler &cgh){
        auto ssha_t_accessor = ssha_t_buffer.get_access<access::mode::write>(cgh);

        auto sshn_u_accessor = sshn_u_buffer.get_access<access::mode::read>(cgh);
        auto sshn_v_accessor = sshn_v_buffer.get_access<access::mode::read>(cgh);
        auto sshn_t_accessor = sshn_t_buffer.get_access<access::mode::read_write>(cgh);
        auto hu_accessor = hu_buffer.get_access<access::mode::read>(cgh);
        auto hv_accessor = hv_buffer.get_access<access::mode::read>(cgh);
        auto un_accessor = un_buffer.get_access<access::mode::read>(cgh);
        auto vn_accessor = vn_buffer.get_access<access::mode::read>(cgh);

        auto area_t_accessor = area_t_buffer.get_access<access::mode::read>(cgh);

        // Continuity kernel
        cgh.parallel_for(range<2>(internal_ystop, internal_xstop), [=](id<2> idx){
            double rtmp1, rtmp2, rtmp3, rtmp4;
            auto jj = idx[0];
            auto ji = idx[1];
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
            
    // Momentum_u kernel (internal domain u points)
    myqueue.submit([&](handler &cgh){
        auto ua_accessor = ua_buffer.get_access<access::mode::write>(cgh);

        auto sshn_u_accessor = sshn_u_buffer.get_access<access::mode::read>(cgh);
        auto sshn_v_accessor = sshn_v_buffer.get_access<access::mode::read>(cgh);
        auto sshn_t_accessor = sshn_t_buffer.get_access<access::mode::read>(cgh);
        auto hu_accessor = hu_buffer.get_access<access::mode::read>(cgh);
        auto hv_accessor = hv_buffer.get_access<access::mode::read>(cgh);
        auto un_accessor = un_buffer.get_access<access::mode::read>(cgh);
        auto vn_accessor = vn_buffer.get_access<access::mode::read>(cgh);
        auto ht_accessor = ht_buffer.get_access<access::mode::read>(cgh);
        auto ssha_u_accessor = ssha_u_buffer.get_access<access::mode::read>(cgh);

        auto tmask_accessor = tmask_buffer.get_access<access::mode::read>(cgh);
        auto area_u_accessor = area_u_buffer.get_access<access::mode::read>(cgh);
        auto dx_u_accessor = dx_u_buffer.get_access<access::mode::read>(cgh);
        auto dx_v_accessor = dx_v_buffer.get_access<access::mode::read>(cgh);
        auto dx_t_accessor = dx_t_buffer.get_access<access::mode::read>(cgh);
        auto dy_u_accessor = dy_u_buffer.get_access<access::mode::read>(cgh);
        auto dy_t_accessor = dy_t_buffer.get_access<access::mode::read>(cgh);
        auto gphiu_accessor = gphiu_buffer.get_access<access::mode::read>(cgh);

        cgh.parallel_for(range<2>(internal_ystop, internal_xstop - 1), [=](id<2> idx){
            double u_e, u_w, v_n, v_s;
            double v_nc, v_sc;
            double depe, depw, deps, depn;
            double hpg, adv, cor, vis;
            double dudx_e, dudx_w;
            double dudy_s, dudy_n;
            double uu_e, uu_n, uu_s, uu_w;

            auto ji = idx[0];
            auto jj = idx[1];
            if (ji < internal_xstart) return;
            if (jj < internal_ystart) return;

            if(tmask_accessor[{jj, ji}] + tmask_accessor[{jj, ji + 1}] <= 0){
                return; // jump over non-computational domain
            }
            if(tmask_accessor[{jj, ji}] <= 0 || tmask_accessor[{jj, ji + 1}] <= 0){
                return; // jump over boundary u
            }

            u_e  = 0.5 * (un_accessor[{jj, ji}] + un_accessor[{jj, ji + 1}]) *
                dy_t_accessor[{jj, ji + 1}];   // add length scale.
            depe = ht_accessor[{jj, ji + 1}] + sshn_t_accessor[{jj, ji + 1}];

            u_w  = 0.5 * (un_accessor[{jj, ji}] + un_accessor[{jj, ji - 1}]) *
                dy_t_accessor[{jj, ji}];   // add length scale
            depw = ht_accessor[{jj, ji}] + sshn_t_accessor[{jj, ji}];

            v_sc = 0.5 * (vn_accessor[{jj - 1, ji}] + vn_accessor[{jj - 1, ji + 1}]);
            v_s  = 0.5 * v_sc * (dx_v_accessor[{jj - 1, ji}] + dx_v_accessor[{jj - 1, ji + 1}]);
            deps = 0.5 * (hv_accessor[{jj - 1, ji}] + sshn_v_accessor[{jj - 1, ji}] +
                hv_accessor[{jj - 1, ji + 1}] + sshn_v_accessor[{jj - 1, ji + 1}]);

            v_nc = 0.5 * (vn_accessor[{jj, ji}] + vn_accessor[{jj, ji + 1}]);
            v_n  = 0.5 * v_nc * (dx_v_accessor[{jj, ji}] + dx_v_accessor[{jj, ji + 1}]);
            depn = 0.5 * (hv_accessor[{jj, ji}] + sshn_v_accessor[{jj, ji}] +
                hv_accessor[{jj, ji + 1}] + sshn_v_accessor[{jj, ji + 1}]);

            // -advection (currently first order upwind)
            uu_w = (0.5 - cl::sycl::copysign(0.5, u_w)) * un_accessor[{jj, ji}] +
                (0.5 + cl::sycl::copysign(0.5, u_w)) * un_accessor[{jj, ji - 1}] ;
            uu_e = (0.5 + cl::sycl::copysign(0.5, u_e)) * un_accessor[{jj, ji}] +
                (0.5 - cl::sycl::copysign(0.5, u_e)) * un_accessor[{jj, ji + 1}] ;

            if(tmask_accessor[{jj - 1, ji}] <=0 || tmask_accessor[{jj - 1, ji + 1}] <= 0){
                uu_s = (0.5 - cl::sycl::copysign(0.5, v_s)) * un_accessor[{jj, ji}];
            }else{
                uu_s = (0.5 - cl::sycl::copysign(0.5, v_s)) * un_accessor[{jj, ji}] +
                    (0.5 + cl::sycl::copysign(0.5, v_s)) * un_accessor[{jj - 1, ji}] ;
            }

            if(tmask_accessor[{jj + 1, ji}] <=0 || tmask_accessor[{jj + 1 + 1, ji}] <= 0){
                uu_n = (0.5 + cl::sycl::copysign(0.5, v_n)) * un_accessor[{jj, ji}];
            }else{
                uu_n = (0.5 + cl::sycl::copysign(0.5, v_n)) * un_accessor[{jj, ji}] +
                    (0.5 - cl::sycl::copysign(0.5, v_n)) * un_accessor[{jj + 1, ji}];
            }

            adv = uu_w * u_w * depw - uu_e * u_e * depe +
                uu_s * v_s * deps - uu_n * v_n * depn;

            // -viscosity

            dudx_e = (un_accessor[{jj, ji + 1}] - un_accessor[{jj, ji}]) / dx_t_accessor[{jj, ji + 1}] *
                (ht_accessor[{jj, ji + 1}] + sshn_t_accessor[{jj, ji + 1}]);
            dudx_w = (un_accessor[{jj, ji}] - un_accessor[{jj, ji - 1}]) / dx_t_accessor[{jj, ji}] *
                (ht_accessor[{jj, ji}] + sshn_t_accessor[{jj, ji}]);
            if(tmask_accessor[{jj - 1, ji}] <=0 || tmask_accessor[{jj - 1, ji + 1}] <= 0){
                dudy_s = 0.0; // slip boundary
            }else{
                dudy_s = (un_accessor[{jj, ji}] - un_accessor[{jj - 1, ji}]) / (dy_u_accessor[{jj, ji}] +
                    dy_u_accessor[{jj - 1, ji}]) * (hu_accessor[{jj, ji}] + sshn_u_accessor[{jj, ji}] +
                    hu_accessor[{jj - 1, ji}] + sshn_u_accessor[{jj - 1, ji}]);
            }

            if(tmask_accessor[{jj + 1, ji}] <= 0 || tmask_accessor[{jj + 1 + 1, ji}] <= 0){
                dudy_n = 0.0; // slip boundary
            }else{
                dudy_n = (un_accessor[{jj + 1, ji}] - un_accessor[{jj, ji}]) / (dy_u_accessor[{jj, ji}] +
                        dy_u_accessor[{jj + 1, ji}]) * (hu_accessor[{jj, ji}] + sshn_u_accessor[{jj, ji}] +
                        hu_accessor[{jj + 1, ji}] + sshn_u_accessor[{jj + 1, ji}]);
            }

            vis = (dudx_e - dudx_w ) * dy_u_accessor[{jj, ji}]  + 
                (dudy_n - dudy_s ) * dx_u_accessor[{jj, ji}] * 0.5;
            vis = visc * vis;   // visc will be an array visc(1:jpijglou) 
                                // for variable viscosity, such as turbulent viscosity

            // -Coriolis' force (can be implemented implicitly)
            cor = 0.5 * (2. * omega * sin(gphiu_accessor[{jj, ji}] * d2r) * (v_sc + v_nc)) * 
                area_u_accessor[{jj, ji}] * (hu_accessor[{jj, ji}] + sshn_u_accessor[{jj, ji}]);

            // -pressure gradient
            hpg = -g * (hu_accessor[{jj, ji}] + sshn_u_accessor[{jj, ji}]) * dy_u_accessor[{jj, ji}] * 
                (sshn_t_accessor[{jj, ji + 1}] - sshn_t_accessor[{jj, ji}]);

            // -linear bottom friction (implemented implicitly.
            ua_accessor[{jj, ji}] = (un_accessor[{jj, ji}] * (hu_accessor[{jj, ji}] + sshn_u_accessor[{jj, ji}]) +
                rdt * (adv + vis + cor + hpg) / area_u_accessor[{jj, ji}]) / (hu_accessor[{jj, ji}] +
                ssha_u_accessor[{jj, ji}]) / (1.0 + cbfr * rdt) ;

        });
    });

    // Momentum_v kernel (internal domain v points)
    for(int jj = internal_ystart; jj <= internal_ystop - 1; jj++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            momentum_v_code(ji, jj, width, va, un, vn, hu, hv, ht, ssha_v, \
                sshn_t, sshn_u, sshn_v, tmask, dx_v, dx_t, dy_u, dy_v, dy_t, \
                area_v, gphiv, rdt, cbfr, visc, omega, d2r, g);
        }
    }

    // Momentum_v kernel (internal domain v points)
    myqueue.submit([&](handler &cgh){
        auto va_accessor = va_buffer.get_access<access::mode::read_write>(cgh);

        auto sshn_u_accessor = sshn_u_buffer.get_access<access::mode::read_write>(cgh);
        auto sshn_v_accessor = sshn_v_buffer.get_access<access::mode::read_write>(cgh);
        auto sshn_t_accessor = sshn_t_buffer.get_access<access::mode::read_write>(cgh);
        auto hu_accessor = hu_buffer.get_access<access::mode::read>(cgh);
        auto hv_accessor = hv_buffer.get_access<access::mode::read>(cgh);
        auto un_accessor = un_buffer.get_access<access::mode::read_write>(cgh);
        auto vn_accessor = vn_buffer.get_access<access::mode::read_write>(cgh);
        auto ht_accessor = ht_buffer.get_access<access::mode::read>(cgh);
        auto ssha_v_accessor = ssha_v_buffer.get_access<access::mode::read>(cgh);

        auto tmask_accessor = tmask_buffer.get_access<access::mode::read>(cgh);
        auto area_v_accessor = area_v_buffer.get_access<access::mode::read>(cgh);
        auto dx_v_accessor = dx_v_buffer.get_access<access::mode::read>(cgh);
        auto dx_t_accessor = dx_t_buffer.get_access<access::mode::read>(cgh);
        auto dy_u_accessor = dy_u_buffer.get_access<access::mode::read>(cgh);
        auto dy_v_accessor = dy_v_buffer.get_access<access::mode::read>(cgh);
        auto dy_t_accessor = dy_t_buffer.get_access<access::mode::read>(cgh);
        auto gphiv_accessor = gphiv_buffer.get_access<access::mode::read>(cgh);

        cgh.parallel_for(range<2>(internal_ystop - 1, internal_xstop), [=](id<2> idx){

            double u_e, u_w, v_n, v_s;
            double u_ec, u_wc, vv_e, vv_n, vv_s, vv_w;
            double depe, depw, deps, depn;
            double hpg, adv, cor, vis;
            double dvdx_e, dvdx_w, dvdy_n, dvdy_s;
            auto ji = idx[0];
            auto jj = idx[1];
            if (ji < internal_xstart) return;
            if (jj < internal_ystart) return;

        });
    });
    myqueue.wait();

#ifdef USE_TIMER
    TimerStop();
    TimerReport();
#endif


    // Boundary conditions bc_ssh kernel (internal domain)
    myqueue.submit([&](handler &cgh){
        auto ssha_t_accessor = ssha_t_buffer.get_access<access::mode::read_write>(cgh);
        auto tmask_accessor = tmask_buffer.get_access<access::mode::read>(cgh);

        cgh.parallel_for(range<2>(internal_ystop, internal_xstop), [=](id<2> idx){
            auto jj = idx[0];
            auto ji = idx[1];
            if (ji < internal_xstart) return;
            if (jj < internal_ystart) return;

            double amp_tide, omega_tide, rtime;

            amp_tide = 0.2;
            omega_tide = 2.0 * 3.14159 / (12.42 * 3600.0);
            rtime = istep * rdt;
  
            if(tmask_accessor[{jj, ji}] <= 0) return;
  
            if(tmask_accessor[{jj-1, ji}] < 0){
                ssha_t_accessor[{jj, ji}] = amp_tide * sin(omega_tide * rtime);
            }else if(tmask_accessor[{jj+1, ji}] < 0){
                ssha_t_accessor[{jj, ji}] = amp_tide * sin(omega_tide * rtime);
            }else if(tmask_accessor[{jj, ji+1}] < 0){
                ssha_t_accessor[{jj, ji}] = amp_tide * sin(omega_tide * rtime);
            }else if(tmask_accessor[{jj, ji-1}] < 0){
                ssha_t_accessor[{jj, ji}] = amp_tide * sin(omega_tide * rtime);
            }
        });
    });

    // Boundary conditions bc_solid_u kernel (whole domain but top x boundary)
    myqueue.submit([&](handler &cgh){
        auto ua_accessor = ua_buffer.get_access<access::mode::read_write>(cgh);
        auto tmask_accessor = tmask_buffer.get_access<access::mode::read>(cgh);

        cgh.parallel_for(range<2>(internal_ystop + 1, internal_xstop), [=](id<2> idx){
            double rtmp1;

            auto jj = idx[0];
            auto ji = idx[1];
            if (ji < internal_xstart - 1) return;
            if (jj < internal_ystart - 1) return;

            if(tmask_accessor[{jj, ji}] * tmask_accessor[{jj, ji+1}] == 0){
                ua_accessor[{jj, ji}] = 0.0;
            }
        });
    });


    // Boundary conditions bc_solid_v kernel (whole domain but top y boundary)
    myqueue.submit([&](handler &cgh){
        auto va_accessor = va_buffer.get_access<access::mode::read_write>(cgh);
        auto tmask_accessor = tmask_buffer.get_access<access::mode::read>(cgh);

        cgh.parallel_for(range<2>(internal_ystop, internal_xstop + 1), [=](id<2> idx){
            double rtmp1;

            auto jj = idx[0];
            auto ji = idx[1];
            if (ji < internal_xstart - 1) return;
            if (jj < internal_ystart - 1) return;

            if(tmask_accessor[{jj, ji}] * tmask_accessor[{jj+1, ji}] == 0){
                va_accessor[{jj, ji}] = 0.0;
            }
        });
    });

    // Boundary conditions bc_flather_u kernel (whole domain but top x boundary)
    myqueue.submit([&](handler &cgh){
        auto ua_accessor = ua_buffer.get_access<access::mode::read_write>(cgh);
        auto sshn_u_accessor = sshn_u_buffer.get_access<access::mode::read>(cgh);
        auto hu_accessor = hu_buffer.get_access<access::mode::read>(cgh);
        auto tmask_accessor = tmask_buffer.get_access<access::mode::read>(cgh);

        cgh.parallel_for(range<2>(internal_ystop + 1, internal_xstop), [=](id<2> idx){
            double rtmp1;

            auto jj = idx[0];
            auto ji = idx[1];
            if (ji < internal_xstart - 1) return;
            if (jj < internal_ystart - 1) return;

            if(tmask_accessor[{jj, ji}] + tmask_accessor[{jj, ji+1}] <= -1) return;

            if(tmask_accessor[{jj, ji}] < 0){
                ua_accessor[{jj, ji}] = ua_accessor[{jj, ji+1}] + sqrt(g/hu_accessor[{jj, ji}]) *
                    (sshn_u_accessor[{jj, ji}] - sshn_u_accessor[{jj, ji+1}]);
            }else if(tmask_accessor[{jj, ji+1}] < 0){
                ua_accessor[{jj, ji}] = ua_accessor[{jj, ji-1}] + sqrt(g/hu_accessor[{jj, ji}]) *
                    (sshn_u_accessor[{jj, ji}] - sshn_u_accessor[{jj, ji-1}]);
            }
        });
    });

    // Boundary conditions bc_flather_v kernel (whole domain but top y boundary)
    myqueue.submit([&](handler &cgh){
        auto va_accessor = va_buffer.get_access<access::mode::read_write>(cgh);
        auto sshn_v_accessor = sshn_v_buffer.get_access<access::mode::read>(cgh);
        auto hv_accessor = hv_buffer.get_access<access::mode::read>(cgh);
        auto tmask_accessor = tmask_buffer.get_access<access::mode::read>(cgh);

        cgh.parallel_for(range<2>(internal_ystop, internal_xstop + 1), [=](id<2> idx){
            double rtmp1;

            auto jj = idx[0];
            auto ji = idx[1];
            if (ji < internal_xstart - 1) return;
            if (jj < internal_ystart - 1) return;

            if(tmask_accessor[{jj, ji}] + tmask_accessor[{jj+1, ji}] <= -1) return;
    
            if(tmask_accessor[{jj, ji}] < 0){
                va_accessor[{jj, ji}] = va_accessor[{jj+1, ji}] + sqrt(g/hv_accessor[{jj, ji}]) *
                    (sshn_v_accessor[{jj, ji}] - sshn_v_accessor[{jj+1, ji}]);
            }else if(tmask_accessor[{jj+1, ji}] < 0){
                va_accessor[{jj, ji}] = va_accessor[{jj-1, ji}] + sqrt(g/hv_accessor[{jj, ji}]) *
                    (sshn_v_accessor[{jj, ji}] - sshn_v_accessor[{jj-1, ji}]);
            }
        });
    });

    // Copy 'next' fields to 'current' fields (whole domain)
    myqueue.submit([&](handler &cgh){
        auto un_accessor = un_buffer.get_access<access::mode::write>(cgh);
        auto ua_accessor = ua_buffer.get_access<access::mode::read>(cgh);
        cgh.copy(ua_accessor, un_accessor);
    });
    myqueue.submit([&](handler &cgh){
        auto vn_accessor = vn_buffer.get_access<access::mode::write>(cgh);
        auto va_accessor = va_buffer.get_access<access::mode::read>(cgh);
        cgh.copy(va_accessor, vn_accessor);
    });
    myqueue.submit([&](handler &cgh){
        auto sshn_t_accessor = sshn_t_buffer.get_access<access::mode::write>(cgh);
        auto ssha_t_accessor = ssha_t_buffer.get_access<access::mode::read>(cgh);
        cgh.copy(ssha_t_accessor, sshn_t_accessor);
    });

    // Time update kernel (internal domain u points)
    myqueue.submit([&](handler &cgh){
        auto sshn_u_accessor = sshn_u_buffer.get_access<access::mode::read_write>(cgh);
        auto sshn_t_accessor = sshn_t_buffer.get_access<access::mode::read>(cgh);
        auto tmask_accessor = tmask_buffer.get_access<access::mode::read>(cgh);
        auto area_t_accessor = area_t_buffer.get_access<access::mode::read>(cgh);
        auto area_u_accessor = area_u_buffer.get_access<access::mode::read>(cgh);

        cgh.parallel_for(range<2>(internal_ystop, internal_xstop - 1), [=](id<2> idx){
            double rtmp1;

            auto jj = idx[0];
            auto ji = idx[1];
            if (ji < internal_xstart) return;
            if (jj < internal_ystart) return;

            if(tmask_accessor[{jj, ji}] + tmask_accessor[{jj, ji + 1}] <= 0){
                return; // jump over non-computational domain
            }

            if(tmask_accessor[{jj, ji}] * tmask_accessor[{jj, ji + 1}] > 0){
                rtmp1 = area_t_accessor[{jj, ji}] * sshn_t_accessor[{jj, ji}] +
                    area_t_accessor[{jj, ji + 1}] * sshn_t_accessor[{jj, ji + 1}];
                sshn_u_accessor[{jj, ji}] = 0.5 * rtmp1 / area_u_accessor[{jj, ji}];
            }else if(tmask_accessor[{jj, ji}] <= 0){
                sshn_u_accessor[{jj, ji}] = sshn_t_accessor[{jj, ji + 1}];
            }else if(tmask_accessor[{jj, ji + 1}] <= 0){
                sshn_u_accessor[{jj, ji}] = sshn_t_accessor[{jj, ji}];
            }
        });
    });

    // Time update kernel (internal domain v points)
    myqueue.submit([&](handler &cgh){
        auto sshn_v_accessor = sshn_v_buffer.get_access<access::mode::read_write>(cgh);
        auto sshn_t_accessor = sshn_t_buffer.get_access<access::mode::read>(cgh);
        auto tmask_accessor = tmask_buffer.get_access<access::mode::read>(cgh);
        auto area_t_accessor = area_t_buffer.get_access<access::mode::read>(cgh);
        auto area_v_accessor = area_v_buffer.get_access<access::mode::read>(cgh);

        cgh.parallel_for(range<2>(internal_ystop - 1, internal_xstop), [=](id<2> idx){
            double rtmp1;
            auto jj = idx[0];
            auto ji = idx[1];
            if (ji < internal_xstart) return;
            if (jj < internal_ystart) return;

            if((tmask_accessor[{jj, ji}] + tmask_accessor[{jj + 1, ji}]) <= 0){
                return; //jump over non-computational domain
            }
            if((tmask_accessor[{jj, ji}] * tmask_accessor[{jj + 1, ji}]) > 0){
                rtmp1 = area_t_accessor[{jj, ji}] * sshn_t_accessor[{jj, ji}] +
                    area_t_accessor[{jj + 1, ji}] * sshn_t_accessor[{jj + 1, ji}];
                sshn_v_accessor[{jj, ji}] = 0.5 * rtmp1 / area_v_accessor[{jj, ji}];
            }else if(tmask_accessor[{jj, ji}] <= 0){
                sshn_v_accessor[{jj, ji}] = sshn_t_accessor[{jj + 1, ji}];
            }else if(tmask_accessor[{jj + 1, ji}] <= 0){
                sshn_v_accessor[{jj, ji}] = sshn_t_accessor[{jj, ji}];
            }
        });
    });
    myqueue.wait();
}

extern "C" void kokkos_read_from_device(void* from, double* to,
                                        int nx, int ny, int width){

}
