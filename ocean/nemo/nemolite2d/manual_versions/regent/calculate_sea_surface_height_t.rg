import "regent"


require("initialise_grid_points")
require("model_init")

local c = regentlib.c
--This is the FIRST loop.
-- Writes to 2 to N, 2 to M
-- Reads from 1 to N, 1 to M

task calculate_sea_surface_t(sea_surface_after : region(ispace(int2d), uvt_field),
                          sea_surface_now : region(ispace(int2d), uvt_field), 
                          sea_bed_to_mean_sea_level: region(ispace(int2d), uvt_field),
                          velocity_now : region(ispace(int2d), uv_field),
                          grid_region : region(ispace(int2d), grid_fields),
                          rdt : double)
     where sea_surface_after * sea_surface_now,
           sea_surface_after * sea_bed_to_mean_sea_level,
           sea_surface_now * sea_bed_to_mean_sea_level,
     writes( sea_surface_after.t ), reads( sea_surface_now.{u,v,t},
                   sea_bed_to_mean_sea_level.{u,v}, velocity_now.{u,v},
                   grid_region.area_t) do


  for point in sea_surface_after do
    var rtmp1 = (sea_surface_now[point].u + sea_bed_to_mean_sea_level[point].u) 
              * velocity_now[point].u

    var rtmp2 = (sea_surface_now[point+{-1,0}].u 
                 + sea_bed_to_mean_sea_level[point + {-1,0}].u )
              * velocity_now[point + {-1,0}].u

    var rtmp3 = (sea_surface_now[point].v + sea_bed_to_mean_sea_level[point].v)
              * velocity_now[point].v

    var rtmp4 = (sea_surface_now[point + {0,-1}].v 
                + sea_bed_to_mean_sea_level[point + {0,-1}].v)
              * velocity_now[point + {0,-1}].v

    sea_surface_after[point].t = sea_surface_now[point].t + (rtmp2 - rtmp1 + rtmp4 - rtmp3) 
                             * rdt / grid_region[point].area_t

  end
end
