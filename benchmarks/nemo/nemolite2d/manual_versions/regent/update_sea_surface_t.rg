import "regent"


require("model_init")
require("initialise_grid_points")

local c = regentlib.c

sin = regentlib.sin(double)

--This is the FOURTH loop
__demand(__leaf)
task update_sea_surface_t(sea_surface : region(ispace(int2d),uvt_time_field),
                          grid_region : region(ispace(int2d), grid_fields),
                          rdt : double, step : int32 )
     where writes(sea_surface.t_after), reads(grid_region.tmask) do


  var fstep : double = double(step)
  var amp_tide : double = 0.2
  var omega_tide : double = 2.0 * 3.14159 / (12.42 * 3600.0)
  var rtime : double = fstep * rdt
  var tide_value : double = amp_tide * sin(omega_tide * rtime)
__demand(__vectorize)
  for point in sea_surface do
      sea_surface[point].t_after = tide_value
  end


end

