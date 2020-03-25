import "regent"

require("initialise_grid_points")
require("model_init")

local c = regentlib.c

sin = c.sin

--This is the SECOND loop.
--Writes to 2 to N, 2 to M-1
-- Reads from 2 to N+1, 2 to M
task update_velocity_ufield(velocity_after: region(ispace(int2d), uv_field),
                            grid_region : region(ispace(int2d), grid_fields),
                            velocity_now : region(ispace(int2d), uv_field),
                            sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
                            sea_surface_now : region(ispace(int2d), uvt_field),
                            sea_surface_after : region(ispace(int2d),uvt_field),
                            visc : double,
                            omega : double,
                            g : double,
                            rdt : double,
                            cbfr : double,
                            d2r : double )
     where writes(velocity_after.u), reads(grid_region.{tmask, dx_t, dy_t, dx_v,
                                                        dy_u, dx_u, gphi_u, area_u},
                                           velocity_now.{u,v},
                       sea_bed_to_mean_sea_level.{u,t,v}, sea_surface_now.{u,v,t},
                       sea_surface_after.u) do
 
   for point in velocity_after do

     if ( grid_region[point].tmask == int1d(1) and grid_region[point + {1,0}].tmask == int1d(1) ) then
-- ! advection
--        u_e  = 0.5 * (un%data(ji,jj) + un%data(ji+1,jj)) * un%grid%dy_t(ji+1,jj)   !add length scale.
--        depe = ht%data(ji+1,jj) + sshn_t%data(ji+1,jj)
        var u_e = 0.5 * (velocity_now[point].u + velocity_now[point + {1,0}].u) 
                * grid_region[point].dy_t
        var depe = sea_bed_to_mean_sea_level[point + {1,0}].t 
                  + sea_surface_now[point + {1,0}].t 

--        u_w  = 0.5 * (un%data(ji,jj) + un%data(ji-1,jj)) * un%grid%dy_t(ji,jj)     !add length scale
--        depw = ht%data(ji,jj) + sshn_t%data(ji,jj)
        var u_w = 0.5 * (velocity_now[point].u + velocity_now[point + {-1,0}].u)
                * grid_region[point].dy_t
        var depw = sea_bed_to_mean_sea_level[point].t
                  + sea_surface_now[point].t
     
--        v_sc = 0.5_go_wp * (vn%data(ji,jj-1) + vn%data(ji+1,jj-1))
--        v_s  = 0.5_go_wp * v_sc * (un%grid%dx_v(ji,jj-1) + un%grid%dx_v(ji+1,jj-1))
--        deps = 0.5_go_wp * (hv%data(ji,jj-1) + sshn_v%data(ji,jj-1) + hv%data(ji+1,jj-1) + &
--               sshn_v%data(ji+1,jj-1))         
        var v_sc = 0.5 * (velocity_now[point + {0,-1}].v + velocity_now[point + {1,-1}].v ) 
        var v_s = 0.5 * v_sc 
                 * ( grid_region[point + {0,-1}].dx_v + grid_region[point + {1, -1}].dx_v )
        var deps = 0.5 * (sea_bed_to_mean_sea_level[point + {0,-1}].v 
                           + sea_surface_now[point + {0,-1}].v
                           + sea_bed_to_mean_sea_level[point + {1,-1}].v 
                           + sea_surface_now[point + {1,-1}].v )

--        v_nc = 0.5_go_wp * (vn%data(ji,jj) + vn%data(ji+1,jj))
--        v_n  = 0.5_go_wp * v_nc * (un%grid%dx_v(ji,jj) + un%grid%dx_v(ji+1,jj))
--        depn = 0.5_go_wp * (hv%data(ji,jj) + sshn_v%data(ji,jj) + hv%data(ji+1,jj) + &
--                     sshn_v%data(ji+1,jj))
        var v_nc = 0.5 * (velocity_now[point].v + velocity_now[point + {1,0}].v)
        var v_n = 0.5 * v_nc * (grid_region[point].dx_v + grid_region[point + {1,0}].dx_v)
        var depn = 0.5 * (sea_bed_to_mean_sea_level[point].v 
                          + sea_surface_now[point].v
                          + sea_bed_to_mean_sea_level[point + {1,0}].v
                          + sea_surface_now[point + {1,0}].v)
     
--    ! -advection (currently first order upwind)
--        uu_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * un%data(ji,jj)              + &
--             & (0.5_go_wp + SIGN(0.5_go_wp, u_w)) * un%data(ji-1,jj)
      var sign_u_w : double = 0.0
      if(u_w >= 0.0) then
        sign_u_w = 1.0 
      else
        sign_u_w = -1.0
      end
      var uu_w = (0.5 - 0.5 * sign_u_w) * velocity_now[point].u 
               + (0.5 + 0.5 * sign_u_w) * velocity_now[point + {-1,0}].u

--        uu_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * un%data(ji,jj)              + &
--             & (0.5_go_wp - SIGN(0.5_go_wp, u_e)) * un%data(ji+1,jj)
      var sign_u_e : double = 0.0
      if(u_e >= 0.0) then
        sign_u_e = 1.0 
      else
        sign_u_e = -1.0
      end
      var uu_e = (0.5 + 0.5 * sign_u_e) * velocity_now[point].u
               + (0.5 - 0.5 * sign_u_e) * velocity_now[point + {1,0}].u


      var uu_s : double = 0.0
      var sign_v_s : double = 0.0
      if(v_s >= 0.0) then
        sign_v_s = 1.0 
      else
        sign_v_s = -1.0
      end
--        IF(un%grid%tmask(ji,jj-1) <=0 .OR. un%grid%tmask(ji+1,jj-1) <= 0) THEN 
--           uu_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * un%data(ji,jj)
--        ELSE
--           uu_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * un%data(ji,jj)              + &
--                & (0.5_go_wp + SIGN(0.5_go_wp, v_s)) * un%data(ji,jj-1)
--        END If
      if (grid_region[point + {0,-1}].tmask <= int1d(0) 
          or grid_region[point + {1,-1}].tmask <= int1d(0) ) then
          uu_s = (0.5 - sign_v_s) * velocity_now[point].u

      else
         uu_s = (0.5 - sign_v_s) * velocity_now[point].u
              + (0.5 + sign_v_s) * velocity_now[point + {0,-1}].u
 
      end


      var uu_n : double = 0.0
      var sign_v_n : double = 0.0
      if(v_n >= 0.0) then
        sign_v_n = 1.0 
      else
        sign_v_n = -1.0
      end
--    IF(un%grid%tmask(ji,jj+1) <=0 .OR. un%grid%tmask(ji+1,jj+1) <= 0) THEN
--       uu_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * un%data(ji,jj)
--    ELSE
--       uu_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * un%data(ji,jj)              + &
--            & (0.5_go_wp - SIGN(0.5_go_wp, v_n)) * un%data(ji,jj+1)
--    END IF
--
      if( grid_region[point + {0,1}].tmask <= int1d(0) 
          or grid_region[point + {1,1}].tmask <= int1d(0)) then
          uu_n = (0.5 + 0.5 * sign_v_n) * velocity_now[point].u 
      else
          uu_n = (0.5 + 0.5 * sign_v_n) * velocity_now[point].u 
               + (0.5 - 0.5 * sign_v_n) * velocity_now[point + {0,1}].u
      end


--    adv = uu_w * u_w * depw - uu_e * u_e * depe + &
--          uu_s * v_s * deps - uu_n * v_n * depn
      var adv = uu_w * u_w * depw - uu_e * u_e * depe 
              + uu_s * v_s * deps - uu_n * v_n * depn

-- ! visocity

--    !kernel  u vis
--    dudx_e = (un%data(ji+1,jj) - un%data(ji,  jj)) / un%grid%dx_t(ji+1,jj) * &
--             (ht%data(ji+1,jj) + sshn_t%data(ji+1,jj))
      regentlib.assert(grid_region[point + {1,0}].dx_t ~= 0.0, "Divide by 0")
      var dudx_e = (velocity_now[point+{1,0}].u - velocity_now[point].u) 
                 / grid_region[point + {1,0}].dx_t 
                 * (sea_bed_to_mean_sea_level[point + {1,0}].t 
                 + sea_surface_now[point + {1,0}].t)
  
--    dudx_w = (un%data(ji,  jj) - un%data(ji-1,jj)) / un%grid%dx_t(ji,  jj) * &
--             (ht%data(ji,  jj) + sshn_t%data(ji,  jj))
      regentlib.assert(grid_region[point].dx_t ~= 0.0, "Divide by 0")
      var dudx_w = (velocity_now[point].u - velocity_now[point + {-1,0}].u)
                 / grid_region[point].dx_t
                 * (sea_bed_to_mean_sea_level[point].t
                 + sea_surface_now[point].t)
 
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
          regentlib.assert(grid_region[point].dy_u + grid_region[point + {0,-1}].dy_u ~= 0.0, "Divide by 0")
          dudy_s = (velocity_now[point].u - velocity_now[point + {0,-1}].u)
                 / (grid_region[point].dy_u + grid_region[point + {0,-1}].dy_u)
                 * (sea_bed_to_mean_sea_level[point].u 
                    + sea_surface_now[point].u
                    + sea_bed_to_mean_sea_level[point + {0,-1}].u
                    + sea_surface_now[point + {0,-1}].u )
      end

     var dudy_n : double = 0.0

--    IF(un%grid%tmask(ji,jj+1) <= 0 .OR. un%grid%tmask(ji+1,jj+1) <= 0) THEN
--       dudy_n = 0.0_go_wp ! slip boundary
--    ELSE
--       dudy_n = (un%data(ji,jj+1) - un%data(ji,jj)) / (un%grid%dy_u(ji,jj) + un%grid%dy_u(ji,jj+1)) * &
--            & (hu%data(ji,jj) + sshn_u%data(ji,jj) + hu%data(ji,jj+1) + sshn_u%data(ji,jj+1))
--    END If

     if( grid_region[point + {0,1}].tmask > int1d(0) and grid_region[point + {1,1}].tmask > int1d(0)) then
          regentlib.assert(grid_region[point].dy_u + grid_region[point + {0,1}].dy_u ~= 0.0, "Divide by 0")
       dudy_n = (velocity_now[point + {0,1}].u - velocity_now[point].u)
               / (grid_region[point].dy_u + grid_region[point + {0,1}].dy_u)
               * ( sea_bed_to_mean_sea_level[point].u
                   + sea_surface_now[point].u
                   + sea_bed_to_mean_sea_level[point + {0,1}].u
                   + sea_surface_now[point + {0,1}].u )
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
              * (sea_bed_to_mean_sea_level[point].u + sea_surface_now[point].u)

--    ! -pressure gradient
--    !start kernel hpg
--    hpg = -g * (hu%data(ji,jj) + sshn_u%data(ji,jj)) * un%grid%dy_u(ji,jj) * &
--           (sshn_t%data(ji+1,jj) - sshn_t%data(ji,jj))
--    !end kernel hpg

     var hpg = -g * (sea_bed_to_mean_sea_level[point].u + sea_surface_now[point].u) 
             * grid_region[point].dy_u
             * (sea_surface_now[point + {1,0}].t - sea_surface_now[point].t)

--    ! -linear bottom friction (implemented implicitly.
--    !kernel ua calculation
--    ua%data(ji,jj) = (un%data(ji,jj) * (hu%data(ji,jj) + sshn_u%data(ji,jj)) + rdt * &
--                 (adv + vis + cor + hpg) / un%grid%area_u(ji,jj)) / &
--                (hu%data(ji,jj) + ssha_u%data(ji,jj)) / (1.0_go_wp + cbfr * rdt)

      regentlib.assert( grid_region[point].area_u ~= 0.0, "Divide by 0")
      regentlib.assert(  ( sea_bed_to_mean_sea_level[point].u
                                 + sea_surface_after[point].u ) ~= 0.0, "Divide by 0")
      regentlib.assert( (1.0 + cbfr * rdt) ~= 0.0, "Divide by 0")
      velocity_after[point].u = ( velocity_now[point].u 
                              * (sea_bed_to_mean_sea_level[point].u
                                 + sea_surface_now[point].u))
                              + rdt * (adv + vis + cor + hpg)
                              / grid_region[point].area_u
                              / ( sea_bed_to_mean_sea_level[point].u
                                 + sea_surface_after[point].u )
                              / (1.0 + cbfr * rdt)

     end --If tmask
   end

end










--This is the THIRD loop.
--Writes to 2 to N-1, 2 to M
-- Reads from TODO
task update_velocity_vfield(velocity_after: region(ispace(int2d), uv_field),
                            grid_region : region(ispace(int2d), grid_fields),
                            velocity_now : region(ispace(int2d), uv_field),
                            sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
                            sea_surface_now : region(ispace(int2d), uvt_field),
                            sea_surface_after : region(ispace(int2d),uvt_field),
                            visc : double,
                            omega : double,
                            g : double,
                            rdt : double,
                            cbfr : double,
                            d2r : double )
     where writes(velocity_after.v), reads(grid_region.{tmask, dx_t, dy_u, dy_t,
                                                        dy_v, dx_v, gphi_v, area_v},
                                           velocity_now.{u,v},
                       sea_bed_to_mean_sea_level.{u,v,t}, sea_surface_now.{u,v,t},
                       sea_surface_after.v,
                       velocity_after.v--TODO REMOVE ASSERT
) do

  for point in velocity_after do

    if( grid_region[point].tmask == int1d(1) and grid_region[point + {0,1}].tmask == int1d(1)) then

--    ! kernel v adv
--    v_n  = 0.5 * (vn%data(ji,jj) + vn%data(ji,jj+1)) * vn%grid%dx_t(ji,jj+1)
--    depn = ht%data(ji,jj+1) + sshn_t%data(ji,jj+1)
      var v_n = 0.5 * ( velocity_now[point].v + velocity_now[point + {0,1}].v) 
              * grid_region[point + {0,1}].dx_t
      var depn = sea_bed_to_mean_sea_level[point + {0,1}].t 
                + sea_surface_now[point + {0,1}].t

--    v_s  = 0.5 * (vn%data(ji,jj) + vn%data(ji,jj-1)) * vn%grid%dx_t(ji,jj) 
--    deps = ht%data(ji,jj) + sshn_t%data(ji,jj)
      var v_s = 0.5 * ( velocity_now[point].v + velocity_now[point + {0,-1}].v)
               * grid_region[point].dx_t
      var deps = sea_bed_to_mean_sea_level[point].t + sea_surface_now[point].t


--    u_wc = 0.5_go_wp * (un%data(ji-1,jj) + un%data(ji-1,jj+1))
--    u_w  = 0.5_go_wp * u_wc * (vn%grid%dy_u(ji-1,jj) + vn%grid%dy_u(ji-1,jj+1))
--    depw = 0.50_go_wp * (hu%data(ji-1,jj) + sshn_u%data(ji-1,jj) + &
--                      hu%data(ji-1,jj+1) + sshn_u%data(ji-1,jj+1))

     var u_wc = 0.5 * (velocity_now[point + {-1,0}].u + velocity_now[point + {-1,1}].u)
     var u_w = 0.5 * u_wc * (grid_region[point + {-1,0}].dy_u 
                            + grid_region[point + {-1,1}].dy_u)
     var depw = 0.5 * (sea_bed_to_mean_sea_level[point + {-1,0}].u 
                       + sea_surface_now[point + {-1,0}].u
                       + sea_bed_to_mean_sea_level[point + {-1,1}].u
                       + sea_surface_now[point + {-1,1}].u)

--    u_ec = 0.5_go_wp * (un%data(ji,jj) + un%data(ji,jj+1))
--    u_e  = 0.5_go_wp * u_ec * (vn%grid%dy_u(ji,jj) + vn%grid%dy_u(ji,jj+1))
--    depe = 0.50_go_wp * (hu%data(ji,jj) + sshn_u%data(ji,jj) + &
--                      hu%data(ji,jj+1) + sshn_u%data(ji,jj+1))
      var u_ec = 0.5 * (velocity_now[point].u + velocity_now[point + {0,1}].u)
      var u_e = 0.5 * u_ec * (grid_region[point].dy_u + grid_region[point + {0,1}].dy_u)
      var depe = 0.5 * (sea_bed_to_mean_sea_level[point].u
                        + sea_surface_now[point].u
                        + sea_bed_to_mean_sea_level[point + {0,1}].u
                        + sea_surface_now[point + {0,1}].u)

--    ! -advection (currently first order upwind)i
--    vv_s = (0.5_go_wp - SIGN(0.5_go_wp, v_s)) * vn%data(ji,jj)     + &
--         & (0.5_go_wp + SIGN(0.5_go_wp, v_s)) * vn%data(ji,jj-1)
      var sign_v_s : double = 0.0
      if(v_s >= 0.0) then
        sign_v_s = 1.0 
      else
        sign_v_s = -1.0
      end
      var vv_s = (0.5 - 0.5*sign_v_s) * velocity_now[point].v 
               + (0.5 + 0.5*sign_v_s) * velocity_now[point + {0,-1}].v

--    vv_n = (0.5_go_wp + SIGN(0.5_go_wp, v_n)) * vn%data(ji,jj)     + &
--         & (0.5_go_wp - SIGN(0.5_go_wp, v_n)) * vn%data(ji,jj+1)
      var sign_v_n : double = 0.0
      if(v_n >= 0.0) then
        sign_v_n = 1.0 
      else
        sign_v_n = -1.0
      end
      var vv_n = (0.5 + 0.5*sign_v_n) * velocity_now[point].v 
               + (0.5 - 0.5*sign_v_n) * velocity_now[point + {0,1}].v

--    IF(vn%grid%tmask(ji-1,jj) <= 0 .OR. vn%grid%tmask(ji-1,jj+1) <= 0) THEN
--       vv_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * vn%data(ji,jj)
--    ELSE
--       vv_w = (0.5_go_wp - SIGN(0.5_go_wp, u_w)) * vn%data(ji,jj)    + &
--            & (0.5_go_wp + SIGN(0.5_go_wp, u_w)) * vn%data(ji-1,jj)
--    END If
      var vv_w : double = 0.0
      var sign_u_w : double = 0.0
      if(u_w >= 0.0) then
        sign_u_w = 1.0 
      else
        sign_u_w = -1.0
      end
      if( grid_region[point + {-1,0}].tmask <= int1d(0) 
       or grid_region[point + {-1,1}].tmask <= int1d(0)) then
        vv_w = (0.5 - 0.5*sign_u_w) * velocity_now[point].u
      else
        vv_w = (0.5 - 0.5*sign_u_w) * velocity_now[point].u 
             + (0.5 + 0.5*sign_u_w) * velocity_now[point + {-1,0}].u 
      end

      var vv_e : double = 0.0
      var sign_u_e : double = 0.0
      if(u_e >= 0.0) then
        sign_u_e = 1.0 
      else
        sign_u_e = -1.0
      end
--    IF(vn%grid%tmask(ji+1,jj) <= 0 .OR. vn%grid%tmask(ji+1,jj+1) <= 0) THEN
--       vv_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * vn%data(ji,jj)
--    ELSE
--       vv_e = (0.5_go_wp + SIGN(0.5_go_wp, u_e)) * vn%data(ji,jj)  + &
--              (0.5_go_wp - SIGN(0.5_go_wp, u_e)) * vn%data(ji+1,jj)
--    END IF
      if( grid_region[point + {1,0}].tmask <= int1d(0)
       or grid_region[point + {1,1}].tmask <= int1d(0) ) then
         vv_e = (0.5 + 0.5 * sign_u_e) * velocity_now[point].v
      else
         vv_e = (0.5 + 0.5 * sign_u_e) * velocity_now[point].v
              + (0.5 - 0.5 * sign_u_e) * velocity_now[point + {1,0}].v
      end

--    adv = vv_w * u_w * depw - vv_e * u_e * depe + &
--          vv_s * v_s * deps - vv_n * v_n * depn
--
      var adv = vv_w * u_w * depw - vv_e * u_e * depe 
              + vv_s * v_s * deps - vv_n * v_n * depn

      if(point == int2d({2,2})) then
        c.printf("-----------ADVECTION----------\n")
        c.printf("vv_w = %19.15e, u_w = %19.15e, depw = %19.15f\n", vv_w, u_w, depw)
        c.printf("vv_e = %19.15e, u_e = %19.15e, depe = %19.15f\n", vv_e, u_e, depe)
        c.printf("vv_s = %19.15e, u_s = %19.15e, deps = %19.15f\n", vv_s, v_s, deps)
        c.printf("vv_n = %19.15e, u_n = %19.15e, depn = %19.15f\n", vv_n, v_n, depn)
        c.printf("--------END ADVECTION---------\n")
      end
--    !end kernel v adv



--    ! -viscosity
--
--
--    !kernel v dis
--    dvdy_n = (vn%data(ji,jj+1) - vn%data(ji,  jj)) / vn%grid%dy_t(ji,jj+1) * &
--                          (ht%data(ji,jj+1) + sshn_t%data(ji,jj+1))
      regentlib.assert( grid_region[point + {0,1}].dy_t ~= 0.0, "Divide by 0")
      var dvdy_n = (velocity_now[point + {0,1}].v - velocity_now[point].v)
                 / grid_region[point + {0,1}].dy_t
                 * (sea_bed_to_mean_sea_level[point].t + sea_surface_now[point].t)

--    dvdy_s = (vn%data(ji,  jj) - vn%data(ji,jj-1)) / vn%grid%dy_t(ji,  jj) * &
--                          (ht%data(ji,  jj) + sshn_t%data(ji,  jj))
      regentlib.assert( grid_region[point].dy_t ~= 0.0, "Divide by 0")
      var dvdy_s = (velocity_now[point].v - velocity_now[point + {0,-1}].v)
                  / grid_region[point].dy_t
                  * (sea_bed_to_mean_sea_level[point].t + sea_surface_now[point].t)


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
         regentlib.assert((grid_region[point].dx_v + grid_region[point + {-1,0}].dx_v) ~= 0.0, "Divide by 0")
         dvdx_w = (velocity_now[point].v - velocity_now[point + {-1,0}].v)
                 / (grid_region[point].dx_v + grid_region[point + {-1,0}].dx_v)
                 * ( sea_bed_to_mean_sea_level[point].v 
                   + sea_surface_now[point].v
                   + sea_bed_to_mean_sea_level[point + {-1,0}].v
                   + sea_surface_now[point + {-1,0}].v)
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
         regentlib.assert((grid_region[point].dx_v + grid_region[point + {1,0}].dx_v) ~= 0.0, "Divide by 0")
         dvdx_e = (velocity_now[point + {1,0}].v - velocity_now[point].v)
                / (grid_region[point].dx_v + grid_region[point + {1,0}].dx_v)
                * (sea_bed_to_mean_sea_level[point].v
                  + sea_surface_now[point].v
                  + sea_bed_to_mean_sea_level[point + {1,0}].v
                  + sea_surface_now[point + {1,0}].v)
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
                                            + sea_surface_now[point].v)
--    !end kernel v cor

--    ! -pressure gradient
--    !kernel v hpg
--    hpg = -g * (hv%data(ji,jj) + sshn_v%data(ji,jj)) * vn%grid%dx_v(ji,jj) * &
--           (sshn_t%data(ji,jj+1) - sshn_t%data(ji,jj))
      var hpg = -g * (sea_bed_to_mean_sea_level[point].v + sea_surface_now[point].v)
              * grid_region[point].dx_v
              * (sea_surface_now[point + {0,1}].t - sea_surface_now[point].t)
--    !kernel v hpg

--    ! -linear bottom friction (implemented implicitly.
--    !kernel ua calculation
--    va%data(ji,jj) = (vn%data(ji,jj) * (hv%data(ji,jj) + sshn_v%data(ji,jj)) + &
--                 rdt * (adv + vis + cor + hpg) / vn%grid%area_v(ji,jj) ) / &
--                 ((hv%data(ji,jj) + ssha_v%data(ji,jj))) / (1.0_go_wp + cbfr * rdt)
      regentlib.assert(grid_region[point].area_v ~= 0.0, "Divide by 0")
      regentlib.assert((sea_bed_to_mean_sea_level[point].v
                                + sea_surface_after[point].v) ~= 0.0, "Divide by 0")
      regentlib.assert((1.0 + cbfr * rdt) ~= 0.0, "Divide by 0")
      regentlib.assert( adv == adv, "adv is NaN")
      regentlib.assert( vis == vis, "vis is NaN")
      regentlib.assert( cor == cor, "cor is NaN")
      regentlib.assert( hpg == hpg, "hpg is NaN")
      regentlib.assert( cbfr == cbfr, "cbfr is NaN")
      regentlib.assert( rdt == rdt, "rdt is NaN")
      regentlib.assert( velocity_now[point].v == velocity_now[point].v, "vel now is NaN")
      regentlib.assert( sea_bed_to_mean_sea_level[point].v == sea_bed_to_mean_sea_level[point].v, "Sea  bed is NaN")
      regentlib.assert( sea_surface_now[point].v == sea_surface_now[point].v, "sea surface is NaN")
      regentlib.assert(rdt == rdt, "rdt is NaN")
      regentlib.assert( grid_region[point].area_v == grid_region[point].area_v, "area_v is NaN")
      regentlib.assert(  sea_surface_after[point].v ==  sea_surface_after[point].v, "sea_sruface_after is NaN")
      
      regentlib.assert(velocity_after[point].v == velocity_after[point].v, "NaN before calculation")
      velocity_after[point].v = (velocity_now[point].v
                                * (sea_bed_to_mean_sea_level[point].v 
                                  + sea_surface_now[point].v)
                                + rdt * (adv + vis + cor + hpg) 
                                / grid_region[point].area_v )
                              / (sea_bed_to_mean_sea_level[point].v 
                                + sea_surface_after[point].v)
                              / (1.0 + cbfr * rdt)
      if(point == int2d({2,2})) then
        var string : int8[5000]
        c.sprintf(&string[0], "NaN in calculation for point %i %i\n, v=%f\n bed=%f\n surface=%f\n rdt=%f\n adv=%f\n vis=%f\n cor=%f\n hpg=%f\n areav=%f\n bed=%f\n suraf=%f\n cbfr=%f\n rdt=%f\n, velo=%f\n bed+surf=%16.7f\n rdt=%16.7f\n comps=%16.7f\n area_v=%16.7f\n bed_after=%16.7f\n sums=%16.7f\n %f\n", point.x, point.y, velocity_now[point].v, sea_bed_to_mean_sea_level[point].v, sea_surface_now[point].v, rdt, adv, vis, cor, hpg, grid_region[point].area_v, sea_bed_to_mean_sea_level[point].v, sea_surface_after[point].v, cbfr, rdt, velocity_now[point].v, (sea_bed_to_mean_sea_level[point].v
                                  + sea_surface_now[point].v),  rdt, (adv + vis + cor + hpg),  grid_region[point].area_v, (sea_bed_to_mean_sea_level[point].v
                                + sea_surface_after[point].v), (1.0 + cbfr * rdt))
--        c.printf("%s", string)
--        c.printf("Surface_now = %17.8f\n", sea_surface_now[point].v)
--        c.printf("sea_bed_level = %17.8f\n",  sea_bed_to_mean_sea_level[point].v)
--        c.printf("sum = %17.8f\n", sea_bed_to_mean_sea_level[point].v + sea_surface_now[point].v)
        c.printf("adv = %19.15f\n vis = %19.15f\n cor = %19.15f\n hpg = %19.15f\n", 
                 adv, vis, cor, hpg)
        regentlib.assert(velocity_after[point].v == velocity_after[point].v, string)
      end
    end
  end

end
