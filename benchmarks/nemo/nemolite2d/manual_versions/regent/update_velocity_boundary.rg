import "regent"

require("initialise_grid_points")
require("model_init")
local c = regentlib.c

--This is the FIFTH loop
__demand(__leaf)
task update_uvel_boundary( velocity : region(ispace(int2d), uv_time_field),
                           grid_region : region(ispace(int2d), grid_fields) )
     where writes(velocity.u_after), reads(grid_region.tmask) do

  __demand(__vectorize)
  for point in velocity do
      velocity[point].u_after = 0.0
  end

end

--This is the SIXTH loop
__demand(__leaf)
task update_vvel_boundary( velocity : region(ispace(int2d), uv_time_field),
                           grid_region : region(ispace(int2d), grid_fields) )
     where writes(velocity.v_after), reads(grid_region.tmask) do

  __demand(__vectorize)
  for point in velocity do
      velocity[point].v_after = 0.0
  end


end
