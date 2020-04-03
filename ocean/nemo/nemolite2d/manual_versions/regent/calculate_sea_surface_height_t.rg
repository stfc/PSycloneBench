import "regent"


require("initialise_grid_points")
require("model_init")

local c = regentlib.c
--This is the FIRST loop.
-- Writes to 2 to N, 2 to M
-- Reads from 1 to N, 1 to M

task calculate_sea_surface_t(sea_surface : region(ispace(int2d), uvt_time_field),
                          sea_bed_to_mean_sea_level: region(ispace(int2d), uvt_field),
                          velocity : region(ispace(int2d), uv_time_field),
                          grid_region : region(ispace(int2d), grid_fields),
                          rdt : double, jpiglo : int, jpjglo : int)
     where -- sea_surface * sea_bed_to_mean_sea_level,
     writes( sea_surface.t_after ), reads( sea_surface.{u_now,v_now,t_now},
                   sea_bed_to_mean_sea_level.{u,v}, velocity.{u_now,v_now},
                   grid_region.area_t) do


  for point in sea_surface do
    var rtmp1 = (sea_surface[point].u_now + sea_bed_to_mean_sea_level[point].u) 
              * velocity[point].u_now

    var rtmp2 = (sea_surface[point+{-1,0}].u_now 
                 + sea_bed_to_mean_sea_level[point + {-1,0}].u )
              * velocity[point + {-1,0}].u_now

    var rtmp3 = (sea_surface[point].v_now + sea_bed_to_mean_sea_level[point].v)
              * velocity[point].v_now

    var rtmp4 = (sea_surface[point + {0,-1}].v_now 
                + sea_bed_to_mean_sea_level[point + {0,-1}].v)
              * velocity[point + {0,-1}].v_now

    sea_surface[point].t_after = sea_surface[point].t_now + (rtmp2 - rtmp1 + rtmp4 - rtmp3) 
                             * rdt / grid_region[point].area_t

  end
end
