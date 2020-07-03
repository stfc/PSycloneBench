import "regent"

require("initialise_grid_points")
require("model_init")

local c = regentlib.c

-- This is the NINTH loop
__demand(__leaf)
task update_velocity_and_t_height( velocity : region(ispace(int2d), uv_time_field),
                                   sea_surface : region(ispace(int2d), uvt_time_field)) where 
                                   writes(velocity.{u_now,v_now}, sea_surface.t_now),
                                   reads(velocity.{u_after,v_after}, sea_surface.t_after) 
do

  __demand(__vectorize)
  for point in velocity do
    velocity[point].u_now = velocity[point].u_after
    velocity[point].v_now = velocity[point].v_after
    sea_surface[point].t_now = sea_surface[point].t_after
  end
end

--This is the TENTH loop
__demand(__leaf)
task update_u_height_launcher( sea_surface : region(ispace(int2d), uvt_time_field),
                               grid_region : region(ispace(int2d), grid_fields), 
                               sea_surface_halos : region(ispace(int2d), uvt_time_field) ) where 
                               writes(sea_surface.u_now), 
                               reads(grid_region.area_u, grid_region.area_t, grid_region.tmask, sea_surface_halos.t_after) 
do

--If statement currently prevents vectorisation
  for point in sea_surface do
    if( grid_region[point].tmask + grid_region[point+{1,0}].tmask > int1d(0)) then
        if(grid_region[point].tmask * grid_region[point + {1,0}].tmask > int1d(0)) then
           var rtmp1 : double = grid_region[point].area_t * sea_surface_halos[point].t_after
                          + grid_region[point+{1,0}].area_t * sea_surface_halos[point+{1,0}].t_after
           sea_surface[point].u_now = 0.5 * rtmp1 / grid_region[point].area_u
        elseif( grid_region[point].tmask <= int1d(0) ) then
           sea_surface[point].u_now = sea_surface_halos[point + {1,0}].t_after
        elseif( grid_region[point + {1,0}].tmask <= int1d(0)) then
           sea_surface[point].u_now = sea_surface_halos[point].t_after
        end
    end
  end


end


--This is the ELEVENTH loop
__demand(__leaf)
task update_v_height_launcher( sea_surface : region(ispace(int2d), uvt_time_field),
                               grid_region : region(ispace(int2d), grid_fields),
                               sea_surface_halos : region(ispace(int2d), uvt_time_field) ) where
                               writes(sea_surface.v_now),
                               reads(grid_region.area_v, grid_region.area_t, grid_region.tmask, sea_surface_halos.t_after) 
do


--If statement currently prevents vectorisation
  for point in sea_surface do
    if( grid_region[point].tmask + grid_region[point+{0,1}].tmask > int1d(0)) then
         if(grid_region[point].tmask * grid_region[point + {0,1}].tmask > int1d(0)) then
           var rtmp1 : double = grid_region[point].area_t * sea_surface_halos[point].t_after
                          + grid_region[point+{0,1}].area_t * sea_surface_halos[point+{0,1}].t_after
           sea_surface[point].v_now = 0.5 * rtmp1 / grid_region[point].area_v
        elseif( grid_region[point].tmask <= int1d(0) ) then
           sea_surface[point].v_now = sea_surface_halos[point + {0,1}].t_after
        elseif( grid_region[point + {0,1}].tmask <= int1d(0)) then
           sea_surface[point].v_now = sea_surface_halos[point].t_after
        end 
    end
  end
end 
