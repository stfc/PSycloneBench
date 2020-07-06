import "regent"

require("initialise_grid_points")
require("model_init")

local c = regentlib.c
local sin = regentlib.sin(double)

-- Regent's inbuilt copysign is currently significantly slower than our own implementation (github.com/StanfordLegion/legion #864)
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

--This is the SECOND loop. This is equivalent to the momentum kernel in the Fortran version.
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
     where 
           writes(velocity.u_after), reads(grid_region.{tmask, dx_t, dy_t, dx_v,
                                                        dy_u, dx_u, gphi_u, area_u},
                                           velocity_full.{u_now,v_now},
                       sea_bed_to_mean_sea_level.{u,t,v}, sea_surface.{u_now,v_now,t_now},
                       sea_surface.u_after) do
--This loop cannot yet vectorize due to copysign implementation (github.com/StanfordLegion/legion #864) and some remaining if statements (#859)
--    __demand(__vectorize)
    for point in velocity do

-- ! advection
        var u_e = 0.5 * (velocity_full[point].u_now + velocity_full[point + {1,0}].u_now) 
                * grid_region[point].dy_t
        var depe = sea_bed_to_mean_sea_level[point + {1,0}].t 
                  + sea_surface[point + {1,0}].t_now 

        var u_w = 0.5 * (velocity_full[point].u_now + velocity_full[point + {-1,0}].u_now)
                * grid_region[point].dy_t
        var depw = sea_bed_to_mean_sea_level[point].t
                  + sea_surface[point].t_now
     
        var v_sc = 0.5 * (velocity_full[point + {0,-1}].v_now + velocity_full[point + {1,-1}].v_now ) 
        var v_s = 0.5 * v_sc 
                 * ( grid_region[point + {0,-1}].dx_v + grid_region[point + {1, -1}].dx_v )
        var deps = 0.5 * (sea_bed_to_mean_sea_level[point + {0,-1}].v 
                           + sea_surface[point + {0,-1}].v_now
                           + sea_bed_to_mean_sea_level[point + {1,-1}].v 
                           + sea_surface[point + {1,-1}].v_now )

        var v_nc = 0.5 * (velocity_full[point].v_now + velocity_full[point + {1,0}].v_now)
        var v_n = 0.5 * v_nc * (grid_region[point].dx_v + grid_region[point + {1,0}].dx_v)
        var depn = 0.5 * (sea_bed_to_mean_sea_level[point].v 
                          + sea_surface[point].v_now
                          + sea_bed_to_mean_sea_level[point + {1,0}].v
                          + sea_surface[point + {1,0}].v_now)
     
--    ! -advection (currently first order upwind)
      var sign_uw = copysign(0.5, u_w)
      var uu_w = (0.5 - sign_uw) * velocity_full[point].u_now 
               + (0.5 + sign_uw) * velocity_full[point + {-1,0}].u_now

      var sign_ue = copysign(0.5, u_e)
      var uu_e = (0.5 + sign_ue) * velocity_full[point].u_now
               + (0.5 - sign_ue) * velocity_full[point + {1,0}].u_now


      var uu_s : double = 0.0
      var sign_vs = copysign(0.5, v_s)
      if (grid_region[point + {0,-1}].tmask <= int1d(0) 
          or grid_region[point + {1,-1}].tmask <= int1d(0) ) then
          uu_s = (0.5 - sign_vs) * velocity_full[point].u_now

      else
         uu_s = (0.5 - sign_vs) * velocity_full[point].u_now
              + (0.5 + sign_vs) * velocity_full[point + {0,-1}].u_now
 
      end


      var uu_n : double = 0.0
      var sign_vn = copysign(0.5, v_n)
      if( grid_region[point + {0,1}].tmask <= int1d(0) 
          or grid_region[point + {1,1}].tmask <= int1d(0)) then
          uu_n = (0.5 + sign_vn) * velocity_full[point].u_now 
      else
          uu_n = (0.5 + sign_vn) * velocity_full[point].u_now 
               + (0.5 - sign_vn) * velocity_full[point + {0,1}].u_now
      end


      var adv = uu_w * u_w * depw - uu_e * u_e * depe 
              + uu_s * v_s * deps - uu_n * v_n * depn

-- ! visocity

      var dudx_e = (velocity_full[point+{1,0}].u_now - velocity_full[point].u_now) 
                 / grid_region[point + {1,0}].dx_t 
                 * (sea_bed_to_mean_sea_level[point + {1,0}].t 
                 + sea_surface[point + {1,0}].t_now)
  
      var dudx_w = (velocity_full[point].u_now - velocity_full[point + {-1,0}].u_now)
                 / grid_region[point].dx_t
                 * (sea_bed_to_mean_sea_level[point].t
                 + sea_surface[point].t_now)
 
      var dudy_s : double = 0.0
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

     if( grid_region[point + {0,1}].tmask > int1d(0) and grid_region[point + {1,1}].tmask > int1d(0)) then
       dudy_n = (velocity_full[point + {0,1}].u_now - velocity_full[point].u_now)
               / (grid_region[point].dy_u + grid_region[point + {0,1}].dy_u)
               * ( sea_bed_to_mean_sea_level[point].u
                   + sea_surface[point].u_now
                   + sea_bed_to_mean_sea_level[point + {0,1}].u
                   + sea_surface[point + {0,1}].u_now )
     end

     var vis =  (dudx_e - dudx_w) * grid_region[point].dy_u 
              + (dudy_n - dudy_s) * grid_region[point].dx_u * 0.5
     vis = visc * vis

--    ! -Coriolis' force (can be implemented implicitly)
      var cor = 0.5 * (2.0 * omega * sin(grid_region[point].gphi_u * d2r) * (v_sc + v_nc))
              * grid_region[point].area_u 
              * (sea_bed_to_mean_sea_level[point].u + sea_surface[point].u_now)

--    ! -pressure gradient
--    !start kernel hpg

     var hpg = -g * (sea_bed_to_mean_sea_level[point].u + sea_surface[point].u_now) 
             * grid_region[point].dy_u
             * (sea_surface[point + {1,0}].t_now - sea_surface[point].t_now)

--    ! -linear bottom friction (implemented implicitly.
      velocity[point].u_after = (velocity_full[point].u_now * 
                        (sea_bed_to_mean_sea_level[point].u + sea_surface[point].u_now) + 
                        rdt * (adv + vis + cor + hpg ) / grid_region[point].area_u) 
                              / ( sea_bed_to_mean_sea_level[point].u
                                 + sea_surface[point].u_after )
                              / (1.0 + cbfr * rdt)     
   end

end

--This is the THIRD loop. This is equivalent to the momentum kernel in the Fortran version.
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
     where 
           writes(velocity.v_after), reads(grid_region.{tmask, dx_t, dy_u, dy_t,
                                                        dy_v, dx_v, gphi_v, area_v},
                                           velocity_full.{u_now,v_now},
                       sea_bed_to_mean_sea_level.{u,v,t}, sea_surface.{u_now,v_now,t_now},
                       sea_surface.v_after
) do

--This loop cannot yet vectorize due to copysign implementation (github.com/StanfordLegion/legion #864) and some remaining if statements (#859)
--    __demand(__vectorize)
  for point in velocity do


--    ! kernel v adv
      var v_n = 0.5 * ( velocity_full[point].v_now + velocity_full[point + {0,1}].v_now) 
              * grid_region[point + {0,1}].dx_t
      var depn = sea_bed_to_mean_sea_level[point + {0,1}].t 
                + sea_surface[point + {0,1}].t_now

      var v_s = 0.5 * ( velocity_full[point].v_now + velocity_full[point + {0,-1}].v_now)
               * grid_region[point].dx_t
      var deps = sea_bed_to_mean_sea_level[point].t + sea_surface[point].t_now


     var u_wc = 0.5 * (velocity_full[point + {-1,0}].u_now + velocity_full[point + {-1,1}].u_now)
     var u_w = 0.5 * u_wc * (grid_region[point + {-1,0}].dy_u 
                            + grid_region[point + {-1,1}].dy_u)
     var depw = 0.5 * (sea_bed_to_mean_sea_level[point + {-1,0}].u 
                       + sea_surface[point + {-1,0}].u_now
                       + sea_bed_to_mean_sea_level[point + {-1,1}].u
                       + sea_surface[point + {-1,1}].u_now)

      var u_ec = 0.5 * (velocity_full[point].u_now + velocity_full[point + {0,1}].u_now)
      var u_e = 0.5 * u_ec * (grid_region[point].dy_u + grid_region[point + {0,1}].dy_u)
      var depe = 0.5 * (sea_bed_to_mean_sea_level[point].u
                        + sea_surface[point].u_now
                        + sea_bed_to_mean_sea_level[point + {0,1}].u
                        + sea_surface[point + {0,1}].u_now)

--    ! -advection (currently first order upwind)i
      var sign_vs = copysign(0.5, v_s)
      var vv_s = (0.5 - sign_vs) * velocity_full[point].v_now
               + (0.5 + sign_vs) * velocity_full[point + {0,-1}].v_now

      var sign_vn = copysign(0.5, v_n)
      var vv_n = (0.5 + sign_vn) * velocity_full[point].v_now 
               + (0.5 - sign_vn) * velocity_full[point + {0,1}].v_now

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
      if( grid_region[point + {1,0}].tmask <= int1d(0)
       or grid_region[point + {1,1}].tmask <= int1d(0) ) then
         vv_e = (0.5 + sign_ue) * velocity_full[point].v_now
      else
         vv_e = (0.5 + sign_ue) * velocity_full[point].v_now
              + (0.5 - sign_ue) * velocity_full[point + {1,0}].v_now
      end

      var adv = vv_w * u_w * depw - vv_e * u_e * depe 
              + vv_s * v_s * deps - vv_n * v_n * depn
--    !end kernel v adv



--    ! -viscosity
      var dvdy_n = (velocity_full[point + {0,1}].v_now - velocity_full[point].v_now)
                 / grid_region[point + {0,1}].dy_t
                 * (sea_bed_to_mean_sea_level[point+{0,1}].t + sea_surface[point+{0,1}].t_now)

      var dvdy_s = (velocity_full[point].v_now - velocity_full[point + {0,-1}].v_now)
                  / grid_region[point].dy_t
                  * (sea_bed_to_mean_sea_level[point].t + sea_surface[point].t_now)


      var dvdx_w : double = 0.0
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
      if(grid_region[point + {1,0}].tmask > int1d(0) 
         and grid_region[point + {1,1}].tmask > int1d(0)) then
         dvdx_e = (velocity_full[point + {1,0}].v_now - velocity_full[point].v_now)
                / (grid_region[point].dx_v + grid_region[point + {1,0}].dx_v)
                * (sea_bed_to_mean_sea_level[point].v
                  + sea_surface[point].v_now
                  + sea_bed_to_mean_sea_level[point + {1,0}].v
                  + sea_surface[point + {1,0}].v_now)
       end


      var vis = (dvdy_n - dvdy_s) * grid_region[point].dx_v
              + (dvdx_e - dvdx_w) * grid_region[point].dy_v * 0.5
      vis = visc * vis

--    !end kernel v dis
--    !kernel v cor
      var cor = -0.5* (2.0 * omega * sin(grid_region[point].gphi_v * d2r) * (u_ec + u_wc))
              * grid_region[point].area_v * (sea_bed_to_mean_sea_level[point].v
                                            + sea_surface[point].v_now)
--    !end kernel v cor

--    ! -pressure gradient
--    !kernel v hpg
      var hpg : double = -g * (sea_bed_to_mean_sea_level[point].v + sea_surface[point].v_now)
              * grid_region[point].dx_v
              * (sea_surface[point + {0,1}].t_now - sea_surface[point].t_now)
--    !kernel v hpg

--    ! -linear bottom friction (implemented implicitly.
--    !kernel ua calculation
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
