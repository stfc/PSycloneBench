import "regent"

require("initialise_grid_points")
require("model_init")

local c = regentlib.c
--local math = terralib.includec("math.h")
local sin = regentlib.sin(double)
--terra sin(f : double)
--  return
--end

-- local copysign = regentlib.copysign(double)

struct convert_struct{
  union{
    d : double;
    l : int64;
  }
}
terra double_as_int64(x : double) : int64
 var c : convert_struct
 c.d = x
 return c.l
end

terra int64_as_double(x: int64) : double
  var c : convert_struct
  c.l = x
  return c.d
end

terra copysign( x : double, y: double) : double
  return int64_as_double(  ( double_as_int64(x) and not(1LL << 63)) ^ (double_as_int64(y) and (1LL << 63) ) )
end

task calculate_halo_size( private_bounds: rect2d) : rect2d
  return rect2d({ private_bounds.lo - {1,1}, private_bounds.hi + {1,1} })
end



--This is the SECOND loop.
--Writes to 2 to N, 2 to M-1
-- Reads from 1 to N+1, 1 to M
__demand(__leaf)
task update_velocity_ufield(velocity: region(ispace(int2d), uv_time_field),
                            grid_region : region(ispace(int2d), grid_fields),
                            velocity_full : region(ispace(int2d), uv_time_field),
                            sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
                            sea_surface : region(ispace(int2d), uvt_time_field),
                            visc : double,
                            omega : double,
                            g : double,
                            rdt : double,
                            cbfr : double,
                            d2r : double)
     where -- velocity <= velocity_full,
           writes(velocity.u_after), reads(grid_region.{tmask, dx_t, dy_t, dx_v,
                                                        dy_u, dx_u, gphi_u, area_u},
                                           velocity_full.{u_now,v_now},
                       sea_bed_to_mean_sea_level.{u,t,v}, sea_surface.{u_now,v_now,t_now},
                       sea_surface.u_after) do
--    __demand(__vectorize)
    for point in velocity do

-- ! advection
--        u_e  = 0.5 * (un%data(ji,jj) + un%data(ji+1,jj)) * un%grid%dy_t(ji+1,jj)   !add length scale.
--        depe = ht%data(ji+1,jj) + sshn_t%data(ji+1,jj)
        var u_e = 0.5 * (velocity_full[point].u_now + velocity_full[point + {1,0}].u_now) 
                * grid_region[point].dy_t
        var depe = sea_bed_to_mean_sea_level[point + {1,0}].t 
                  + sea_surface[point + {1,0}].t_now 

--        u_w  = 0.5 * (un%data(ji,jj) + un%data(ji-1,jj)) * un%grid%dy_t(ji,jj)     !add length scale
--        depw = ht%data(ji,jj) + sshn_t%data(ji,jj)
        var u_w = 0.5 * (velocity_full[point].u_now + velocity_full[point + {-1,0}].u_now)
                * grid_region[point].dy_t
        var depw = sea_bed_to_mean_sea_level[point].t
                  + sea_surface[point].t_now
     
--        v_sc = 0.5_go_wp * (vn%data(ji,jj-1) + vn%data(ji+1,jj-1))
--        v_s  = 0.5_go_wp * v_sc * (un%grid%dx_v(ji,jj-1) + un%grid%dx_v(ji+1,jj-1))
--        deps = 0.5_go_wp * (hv%data(ji,jj-1) + sshn_v%data(ji,jj-1) + hv%data(ji+1,jj-1) + &
--               sshn_v%data(ji+1,jj-1))         
        var v_sc = 0.5 * (velocity_full[point + {0,-1}].v_now + velocity_full[point + {1,-1}].v_now ) 
        var v_s = 0.5 * v_sc 
                 * ( grid_region[point + {0,-1}].dx_v + grid_region[point + {1, -1}].dx_v )
        var deps = 0.5 * (sea_bed_to_mean_sea_level[point + {0,-1}].v 
                           + sea_surface[point + {0,-1}].v_now
                           + sea_bed_to_mean_sea_level[point + {1,-1}].v 
                           + sea_surface[point + {1,-1}].v_now )

--        v_nc = 0.5_go_wp * (vn%data(ji,jj) + vn%data(ji+1,jj))
--        v_n  = 0.5_go_wp * v_nc * (un%grid%dx_v(ji,jj) + un%grid%dx_v(ji+1,jj))
--        depn = 0.5_go_wp * (hv%data(ji,jj) + sshn_v%data(ji,jj) + hv%data(ji+1,jj) + &
--                     sshn_v%data(ji+1,jj))
        var v_nc = 0.5 * (velocity_full[point].v_now + velocity_full[point + {1,0}].v_now)
        var v_n = 0.5 * v_nc * (grid_region[point].dx_v + grid_region[point + {1,0}].dx_v)
        var depn = 0.5 * (sea_bed_to_mean_sea_level[point].v 
                          + sea_surface[point].v_now
                          + sea_bed_to_mean_sea_level[point + {1,0}].v
                          + sea_surface[point + {1,0}].v_now)
     
--    ! -advection (currently first order upwind)
--        uu_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * un%data(ji,jj)              + &
--             & (0.5_go_wp + SIGN(0.5_go_wp, u_w)) * un%data(ji-1,jj)
      var sign_uw = copysign(0.5, u_w)
      var uu_w = (0.5 - sign_uw) * velocity_full[point].u_now 
               + (0.5 + sign_uw) * velocity_full[point + {-1,0}].u_now

--        uu_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * un%data(ji,jj)              + &
--             & (0.5_go_wp - SIGN(0.5_go_wp, u_e)) * un%data(ji+1,jj)
      var sign_ue = copysign(0.5, u_e)
      var uu_e = (0.5 + sign_ue) * velocity_full[point].u_now
               + (0.5 - sign_ue) * velocity_full[point + {1,0}].u_now


      var uu_s : double = 0.0
--        IF(un%grid%tmask(ji,jj-1) <=0 .OR. un%grid%tmask(ji+1,jj-1) <= 0) THEN 
--           uu_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * un%data(ji,jj)
--        ELSE
--           uu_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * un%data(ji,jj)              + &
--                & (0.5_go_wp + SIGN(0.5_go_wp, v_s)) * un%data(ji,jj-1)
--        END If
      var sign_vs = copysign(0.5, v_s)
      if (grid_region[point + {0,-1}].tmask <= int1d(0) 
          or grid_region[point + {1,-1}].tmask <= int1d(0) ) then
          uu_s = (0.5 - sign_vs) * velocity_full[point].u_now

      else
         uu_s = (0.5 - sign_vs) * velocity_full[point].u_now
              + (0.5 + sign_vs) * velocity_full[point + {0,-1}].u_now
 
      end


      var uu_n : double = 0.0
--    IF(un%grid%tmask(ji,jj+1) <=0 .OR. un%grid%tmask(ji+1,jj+1) <= 0) THEN
--       uu_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * un%data(ji,jj)
--    ELSE
--       uu_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * un%data(ji,jj)              + &
--            & (0.5_go_wp - SIGN(0.5_go_wp, v_n)) * un%data(ji,jj+1)
--    END IF
--
      var sign_vn = copysign(0.5, v_n)
      if( grid_region[point + {0,1}].tmask <= int1d(0) 
          or grid_region[point + {1,1}].tmask <= int1d(0)) then
          uu_n = (0.5 + sign_vn) * velocity_full[point].u_now 
      else
          uu_n = (0.5 + sign_vn) * velocity_full[point].u_now 
               + (0.5 - sign_vn) * velocity_full[point + {0,1}].u_now
      end


--    adv = uu_w * u_w * depw - uu_e * u_e * depe + &
--          uu_s * v_s * deps - uu_n * v_n * depn
      var adv = uu_w * u_w * depw - uu_e * u_e * depe 
              + uu_s * v_s * deps - uu_n * v_n * depn

-- ! visocity

--    !kernel  u vis
--    dudx_e = (un%data(ji+1,jj) - un%data(ji,  jj)) / un%grid%dx_t(ji+1,jj) * &
--             (ht%data(ji+1,jj) + sshn_t%data(ji+1,jj))
      var dudx_e = (velocity_full[point+{1,0}].u_now - velocity_full[point].u_now) 
                 / grid_region[point + {1,0}].dx_t 
                 * (sea_bed_to_mean_sea_level[point + {1,0}].t 
                 + sea_surface[point + {1,0}].t_now)
  
--    dudx_w = (un%data(ji,  jj) - un%data(ji-1,jj)) / un%grid%dx_t(ji,  jj) * &
--             (ht%data(ji,  jj) + sshn_t%data(ji,  jj))
      var dudx_w = (velocity_full[point].u_now - velocity_full[point + {-1,0}].u_now)
                 / grid_region[point].dx_t
                 * (sea_bed_to_mean_sea_level[point].t
                 + sea_surface[point].t_now)
 
      var dudy_s : double = 0.0
--    IF(un%grid%tmask(ji,jj-1) <=0 .OR. un%grid%tmask(ji+1,jj-1) <= 0) THEN
-- Reversing logic as already set to 0 before the if stataement
--       dudy_s = 0.0_go_wp !slip boundary
--    ELSE
--       dudy_s = (un%data(ji,jj) - un%data(ji,jj-1)) / (un%grid%dy_u(ji,jj) + un%grid%dy_u(ji,jj-1)) * &
--            & (hu%data(ji,jj) + sshn_u%data(ji,jj) + hu%data(ji,jj-1) + sshn_u%data(ji,jj-1))
--    END IF
      if( grid_region[point + {0,-1}].tmask > int1d(0) 
          and grid_region[point + {1,-1}].tmask > int1d(0)) then 
          dudy_s = (velocity_full[point].u_now - velocity_full[point + {0,-1}].u_now)
                 / (grid_region[point].dy_u + grid_region[point + {0,-1}].dy_u)
                 * (sea_bed_to_mean_sea_level[point].u 
                    + sea_surface[point].u_now
                    + sea_bed_to_mean_sea_level[point + {0,-1}].u
                    + sea_surface[point + {0,-1}].u_now )
      end

     var dudy_n : double = 0.0

--    IF(un%grid%tmask(ji,jj+1) <= 0 .OR. un%grid%tmask(ji+1,jj+1) <= 0) THEN
--       dudy_n = 0.0_go_wp ! slip boundary
--    ELSE
--       dudy_n = (un%data(ji,jj+1) - un%data(ji,jj)) / (un%grid%dy_u(ji,jj) + un%grid%dy_u(ji,jj+1)) * &
--            & (hu%data(ji,jj) + sshn_u%data(ji,jj) + hu%data(ji,jj+1) + sshn_u%data(ji,jj+1))
--    END If

     if( grid_region[point + {0,1}].tmask > int1d(0) and grid_region[point + {1,1}].tmask > int1d(0)) then
       dudy_n = (velocity_full[point + {0,1}].u_now - velocity_full[point].u_now)
               / (grid_region[point].dy_u + grid_region[point + {0,1}].dy_u)
               * ( sea_bed_to_mean_sea_level[point].u
                   + sea_surface[point].u_now
                   + sea_bed_to_mean_sea_level[point + {0,1}].u
                   + sea_surface[point + {0,1}].u_now )
     end

--    vis = (dudx_e - dudx_w ) * un%grid%dy_u(ji,jj)  + &
--         & (dudy_n - dudy_s ) * un%grid%dx_u(ji,jj) * 0.5_go_wp
     var vis =  (dudx_e - dudx_w) * grid_region[point].dy_u 
              + (dudy_n - dudy_s) * grid_region[point].dx_u * 0.5
--    vis = visc * vis  
     vis = visc * vis

--    ! -Coriolis' force (can be implemented implicitly)
--    !kernel cor
--    cor = 0.5_go_wp * (2._go_wp * omega * SIN(un%grid%gphiu(ji,jj) * d2r) * (v_sc + v_nc)) * &
--         & un%grid%area_u(ji,jj) * (hu%data(ji,jj) + sshn_u%data(ji,jj))
      var cor = 0.5 * (2.0 * omega * sin(grid_region[point].gphi_u * d2r) * (v_sc + v_nc))
              * grid_region[point].area_u 
              * (sea_bed_to_mean_sea_level[point].u + sea_surface[point].u_now)

--    ! -pressure gradient
--    !start kernel hpg
--    hpg = -g * (hu%data(ji,jj) + sshn_u%data(ji,jj)) * un%grid%dy_u(ji,jj) * &
--           (sshn_t%data(ji+1,jj) - sshn_t%data(ji,jj))
--    !end kernel hpg

     var hpg = -g * (sea_bed_to_mean_sea_level[point].u + sea_surface[point].u_now) 
             * grid_region[point].dy_u
             * (sea_surface[point + {1,0}].t_now - sea_surface[point].t_now)

--    ! -linear bottom friction (implemented implicitly.
--    !kernel ua calculation
--    ua%data(ji,jj) = (un%data(ji,jj) * (hu%data(ji,jj) + sshn_u%data(ji,jj)) + rdt * &
--                 (adv + vis + cor + hpg) / un%grid%area_u(ji,jj)) / &
--                (hu%data(ji,jj) + ssha_u%data(ji,jj)) / (1.0_go_wp + cbfr * rdt)

      velocity[point].u_after = (velocity_full[point].u_now * 
                        (sea_bed_to_mean_sea_level[point].u + sea_surface[point].u_now) + 
                        rdt * (adv + vis + cor + hpg ) / grid_region[point].area_u) 
                              / ( sea_bed_to_mean_sea_level[point].u
                                 + sea_surface[point].u_after )
                              / (1.0 + cbfr * rdt)     
   end

end


--task update_velocity_ufield_launcher( velocities : region(ispace(int2d, uvt_time_field)),
--                            grid_region : region(ispace(int2d), grid_fields),
--                            velocity_halos : region(ispace(int2d, uvt_time_field)),
--                            sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
--                            sea_surface : region(ispace(int2d), uvt_time_field),
--                            visc : double,
--                            omega : double,
--                            g : double,
--                            rdt : double,
--                            cbfr : double,
--                            d2r : double )
--     where sea_bed_to_mean_sea_level * sea_surface,
----           writes(velocities.u_after),
--           reads(grid_region.{tmask, dx_t, dy_t, dx_v,
--                                                        dy_u, dx_u, gphi_u, area_u},
----                                           velocity_halos.{u_now,v_now},
--                       sea_bed_to_mean_sea_level.{u,t,v}, sea_surface.{u_now,v_now,t_now},
--                       sea_surface.u_after) do
----      __demand(__trace, __index_launch)
----      for point in velocities.colors do
----        update_velocity_ufield(velocities[point],
----                               grid_region,
----                               velocity_halos[point],
----                               sea_bed_to_mean_sea_level,
----                               sea_surface,
----                               visc, omega, g, rdt, cbfr, d2r)
----      end
--end








--This is the THIRD loop.
--Writes to 2 to N-1, 2 to M
-- Reads from TODO
task update_velocity_vfield(velocity: region(ispace(int2d), uv_time_field),
                            grid_region : region(ispace(int2d), grid_fields),
                            velocity_full : region(ispace(int2d), uv_time_field),
                            sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
                            sea_surface : region(ispace(int2d), uvt_time_field),
                            visc : double,
                            omega : double,
                            g : double,
                            rdt : double,
                            cbfr : double,
                            d2r : double)
     where --velocity <= velocity_full,
           writes(velocity.v_after), reads(grid_region.{tmask, dx_t, dy_u, dy_t,
                                                        dy_v, dx_v, gphi_v, area_v},
                                           velocity_full.{u_now,v_now},
                       sea_bed_to_mean_sea_level.{u,v,t}, sea_surface.{u_now,v_now,t_now},
                       sea_surface.v_after
) do

  for point in velocity do


--    ! kernel v adv
--    v_n  = 0.5 * (vn%data(ji,jj) + vn%data(ji,jj+1)) * vn%grid%dx_t(ji,jj+1)
--    depn = ht%data(ji,jj+1) + sshn_t%data(ji,jj+1)
      var v_n = 0.5 * ( velocity_full[point].v_now + velocity_full[point + {0,1}].v_now) 
              * grid_region[point + {0,1}].dx_t
      var depn = sea_bed_to_mean_sea_level[point + {0,1}].t 
                + sea_surface[point + {0,1}].t_now

--    v_s  = 0.5 * (vn%data(ji,jj) + vn%data(ji,jj-1)) * vn%grid%dx_t(ji,jj) 
--    deps = ht%data(ji,jj) + sshn_t%data(ji,jj)
      var v_s = 0.5 * ( velocity_full[point].v_now + velocity_full[point + {0,-1}].v_now)
               * grid_region[point].dx_t
      var deps = sea_bed_to_mean_sea_level[point].t + sea_surface[point].t_now


--    u_wc = 0.5_go_wp * (un%data(ji-1,jj) + un%data(ji-1,jj+1))
--    u_w  = 0.5_go_wp * u_wc * (vn%grid%dy_u(ji-1,jj) + vn%grid%dy_u(ji-1,jj+1))
--    depw = 0.50_go_wp * (hu%data(ji-1,jj) + sshn_u%data(ji-1,jj) + &
--                      hu%data(ji-1,jj+1) + sshn_u%data(ji-1,jj+1))

     var u_wc = 0.5 * (velocity_full[point + {-1,0}].u_now + velocity_full[point + {-1,1}].u_now)
     var u_w = 0.5 * u_wc * (grid_region[point + {-1,0}].dy_u 
                            + grid_region[point + {-1,1}].dy_u)
     var depw = 0.5 * (sea_bed_to_mean_sea_level[point + {-1,0}].u 
                       + sea_surface[point + {-1,0}].u_now
                       + sea_bed_to_mean_sea_level[point + {-1,1}].u
                       + sea_surface[point + {-1,1}].u_now)

--    u_ec = 0.5_go_wp * (un%data(ji,jj) + un%data(ji,jj+1))
--    u_e  = 0.5_go_wp * u_ec * (vn%grid%dy_u(ji,jj) + vn%grid%dy_u(ji,jj+1))
--    depe = 0.50_go_wp * (hu%data(ji,jj) + sshn_u%data(ji,jj) + &
--                      hu%data(ji,jj+1) + sshn_u%data(ji,jj+1))
      var u_ec = 0.5 * (velocity_full[point].u_now + velocity_full[point + {0,1}].u_now)
      var u_e = 0.5 * u_ec * (grid_region[point].dy_u + grid_region[point + {0,1}].dy_u)
      var depe = 0.5 * (sea_bed_to_mean_sea_level[point].u
                        + sea_surface[point].u_now
                        + sea_bed_to_mean_sea_level[point + {0,1}].u
                        + sea_surface[point + {0,1}].u_now)

--    ! -advection (currently first order upwind)i
--    vv_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * vn%data(ji,jj)     + &
--         & (0.5_go_wp + SIGN(0.5_go_wp, v_s)) * vn%data(ji,jj-1)
      var sign_vs = copysign(0.5, v_s)
      var vv_s = (0.5 - sign_vs) * velocity_full[point].v_now
               + (0.5 + sign_vs) * velocity_full[point + {0,-1}].v_now

--    vv_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * vn%data(ji,jj)     + &
--         & (0.5_go_wp - SIGN(0.5_go_wp, v_n)) * vn%data(ji,jj+1)
      var sign_vn = copysign(0.5, v_n)
      var vv_n = (0.5 + sign_vn) * velocity_full[point].v_now 
               + (0.5 - sign_vn) * velocity_full[point + {0,1}].v_now

--    IF(vn%grid%tmask(ji-1,jj) <= 0 .OR. vn%grid%tmask(ji-1,jj+1) <= 0) THEN
--       vv_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * vn%data(ji,jj)
--    ELSE
--       vv_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * vn%data(ji,jj)    + &
--            & (0.5_go_wp + SIGN(0.5_go_wp, u_w)) * vn%data(ji-1,jj)
--    END If
      var vv_w : double = 0.0
      var sign_uw = copysign(0.5, u_w)
      if( grid_region[point + {-1,0}].tmask <= int1d(0) 
       or grid_region[point + {-1,1}].tmask <= int1d(0)) then
        vv_w = (0.5 - sign_uw) * velocity_full[point].v_now
      else
        vv_w = (0.5 - sign_uw) * velocity_full[point].v_now
             + (0.5 + sign_uw) * velocity_full[point + {-1,0}].v_now 
      end

      var vv_e : double = 0.0
      var sign_ue = copysign(0.5, u_e)
--    IF(vn%grid%tmask(ji+1,jj) <= 0 .OR. vn%grid%tmask(ji+1,jj+1) <= 0) THEN
--       vv_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * vn%data(ji,jj)
--    ELSE
--       vv_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * vn%data(ji,jj)  + &
--              (0.5_go_wp - SIGN(0.5_go_wp, u_e)) * vn%data(ji+1,jj)
--    END IF
      if( grid_region[point + {1,0}].tmask <= int1d(0)
       or grid_region[point + {1,1}].tmask <= int1d(0) ) then
         vv_e = (0.5 + sign_ue) * velocity_full[point].v_now
      else
         vv_e = (0.5 + sign_ue) * velocity_full[point].v_now
              + (0.5 - sign_ue) * velocity_full[point + {1,0}].v_now
      end

--    adv = vv_w * u_w * depw - vv_e * u_e * depe + &
--          vv_s * v_s * deps - vv_n * v_n * depn
--
      var adv = vv_w * u_w * depw - vv_e * u_e * depe 
              + vv_s * v_s * deps - vv_n * v_n * depn
--    !end kernel v adv



--    ! -viscosity
--
--
--    !kernel v dis
--    dvdy_n = (vn%data(ji,jj+1) - vn%data(ji,  jj)) / vn%grid%dy_t(ji,jj+1) * &
--                          (ht%data(ji,jj+1) + sshn_t%data(ji,jj+1))
      var dvdy_n = (velocity_full[point + {0,1}].v_now - velocity_full[point].v_now)
                 / grid_region[point + {0,1}].dy_t
                 * (sea_bed_to_mean_sea_level[point+{0,1}].t + sea_surface[point+{0,1}].t_now)

--    dvdy_s = (vn%data(ji,  jj) - vn%data(ji,jj-1)) / vn%grid%dy_t(ji,  jj) * &
--                          (ht%data(ji,  jj) + sshn_t%data(ji,  jj))
      var dvdy_s = (velocity_full[point].v_now - velocity_full[point + {0,-1}].v_now)
                  / grid_region[point].dy_t
                  * (sea_bed_to_mean_sea_level[point].t + sea_surface[point].t_now)


      var dvdx_w : double = 0.0
--    IF(vn%grid%tmask(ji-1,jj) <= 0 .OR. vn%grid%tmask(ji-1,jj+1) <= 0) THEN
-- Reversing logic as already set to 0 before the if stataement
--       dvdx_w = 0.0_go_wp !slip boundary
--    ELSE
--       dvdx_w = (vn%data(ji,jj) - vn%data(ji-1,jj)) / &
--                (vn%grid%dx_v(ji,jj) + vn%grid%dx_v(ji-1,jj)) * &
--                (hv%data(ji,jj) + sshn_v%data(ji,jj) + hv%data(ji-1,jj) + sshn_v%data(ji-1,jj))
--    END IF
      if(grid_region[point + {-1,0}].tmask > int1d(0) 
         and grid_region[point + {-1,1}].tmask > int1d(0)) then
         dvdx_w = (velocity_full[point].v_now - velocity_full[point + {-1,0}].v_now)
                 / (grid_region[point].dx_v + grid_region[point + {-1,0}].dx_v)
                 * ( sea_bed_to_mean_sea_level[point].v 
                   + sea_surface[point].v_now
                   + sea_bed_to_mean_sea_level[point + {-1,0}].v
                   + sea_surface[point + {-1,0}].v_now)
      end

      var dvdx_e = 0.0
--    IF(vn%grid%tmask(ji+1,jj) <= 0 .OR. vn%grid%tmask(ji+1,jj+1) <= 0) THEN
--       dvdx_e = 0.0_go_wp ! slip boundary
--    ELSE
--       dvdx_e = (vn%data(ji+1,jj) - vn%data(ji,jj)) / (vn%grid%dx_v(ji,jj) + vn%grid%dx_v(ji+1,jj)) * &
--                  (hv%data(ji,jj) + sshn_v%data(ji,jj) + hv%data(ji+1,jj) + sshn_v%data(ji+1,jj))
--    END If
      if(grid_region[point + {1,0}].tmask > int1d(0) 
         and grid_region[point + {1,1}].tmask > int1d(0)) then
         dvdx_e = (velocity_full[point + {1,0}].v_now - velocity_full[point].v_now)
                / (grid_region[point].dx_v + grid_region[point + {1,0}].dx_v)
                * (sea_bed_to_mean_sea_level[point].v
                  + sea_surface[point].v_now
                  + sea_bed_to_mean_sea_level[point + {1,0}].v
                  + sea_surface[point + {1,0}].v_now)
       end


--    vis = (dvdy_n - dvdy_s ) * vn%grid%dx_v(ji,jj)  + &
--          (dvdx_e - dvdx_w ) * vn%grid%dy_v(ji,jj) * 0.5_go_wp
      var vis = (dvdy_n - dvdy_s) * grid_region[point].dx_v
              + (dvdx_e - dvdx_w) * grid_region[point].dy_v * 0.5
      vis = visc * vis

--    !end kernel v dis
--    !kernel v cor
--    cor = -0.5_go_wp*(2._go_wp * omega * SIN(vn%grid%gphiv(ji,jj) * d2r) * (u_ec + u_wc)) * &
--               vn%grid%area_v(ji,jj) * (hv%data(ji,jj) + sshn_v%data(ji,jj))
      var cor = -0.5* (2.0 * omega * sin(grid_region[point].gphi_v * d2r) * (u_ec + u_wc))
              * grid_region[point].area_v * (sea_bed_to_mean_sea_level[point].v
                                            + sea_surface[point].v_now)
--    !end kernel v cor

--    ! -pressure gradient
--    !kernel v hpg
--    hpg = -g * (hv%data(ji,jj) + sshn_v%data(ji,jj)) * vn%grid%dx_v(ji,jj) * &
--           (sshn_t%data(ji,jj+1) - sshn_t%data(ji,jj))
      var hpg : double = -g * (sea_bed_to_mean_sea_level[point].v + sea_surface[point].v_now)
              * grid_region[point].dx_v
              * (sea_surface[point + {0,1}].t_now - sea_surface[point].t_now)
--    !kernel v hpg

--    ! -linear bottom friction (implemented implicitly.
--    !kernel ua calculation
--    va%data(ji,jj) = (vn%data(ji,jj) * (hv%data(ji,jj) + sshn_v%data(ji,jj)) + &
--                 rdt * (adv + vis + cor + hpg) / vn%grid%area_v(ji,jj) ) / &
--                 ((hv%data(ji,jj) + ssha_v%data(ji,jj))) / (1.0_go_wp + cbfr * rdt)
      velocity[point].v_after = (velocity_full[point].v_now
                                * (sea_bed_to_mean_sea_level[point].v 
                                  + sea_surface[point].v_now)
                                + rdt * (adv + vis + cor + hpg) 
                                / grid_region[point].area_v )
                              / (sea_bed_to_mean_sea_level[point].v 
                                + sea_surface[point].v_after)
                              / (1.0 + cbfr * rdt)
  end

end


--task update_velocity_vfield_launcher( velocity_after: region(ispace(int2d), uv_field),
--                            grid_region : region(ispace(int2d), grid_fields),
--                            velocity_now : region(ispace(int2d), uv_field),
--                            sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
--                            sea_surface_now : region(ispace(int2d), uvt_field),
--                            sea_surface_after : region(ispace(int2d),uvt_field),
--                            visc : double,
--                            omega : double,
--                            g : double,
--                            rdt : double,
--                            cbfr : double,
--                            d2r : double )
--     where velocity_after * velocity_now,
--           sea_bed_to_mean_sea_level * sea_surface_now,
--           sea_bed_to_mean_sea_level * sea_surface_after,
--           sea_surface_now * sea_surface_after,
--           writes(velocity_after.v), reads(grid_region.{tmask, dx_t, dy_u, dy_t,
--                                                        dy_v, dx_v, gphi_v, area_v},
--                                           velocity_now.{u,v},
--                       sea_bed_to_mean_sea_level.{u,v,t}, sea_surface_now.{u,v,t},
--                       sea_surface_after.v
--) do
--
--    var x_dim = velocity_after.bounds.hi.x - velocity_after.bounds.lo.x
--    var y_dim = velocity_after.bounds.hi.y - velocity_after.bounds.lo.y
--    var x_tiles : int32 = x_dim / tilesize
--    var y_tiles : int32 = y_dim / tilesize
--
--    var partition_space = ispace(int2d, {x=x_tiles, y=y_tiles})
--    var partitioned_vel_after = partition(equal, velocity_after, partition_space)

--    var halo_grids = image(grid_region, partitioned_vel_after, calculate_halo_size)
--    var halo_vel_now = image(velocity_now, partitioned_vel_after, calculate_halo_size)
--    var halo_mean_level = image(sea_bed_to_mean_sea_level, partitioned_vel_after, calculate_halo_size)
--    var halo_mean_sn = image(sea_surface_now, partitioned_vel_after, calculate_halo_size)
--    var halo_mean_sa = image(sea_surface_after, partitioned_vel_after, calculate_halo_size)

--    __demand(__index_launch, __trace)
--    for point in partition_space do
--      update_velocity_vfield( partitioned_vel_after[point],
--                              grid_region, velocity_now, sea_bed_to_mean_sea_level, sea_surface_now, sea_surface_after,
--                              --halo_grids[point],
--                              --halo_vel_now[point],
--                              --halo_mean_level[point],
--                              --halo_mean_sn[point],
--                              --halo_mean_sa[point],
--                              visc, omega, g, rdt, cbfr, d2r)
--    end
--    __delete(partitioned_vel_after)

--end
