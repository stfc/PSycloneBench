import "regent"


--requires("model_init")
--requires("initialise_grid_points")

local stdlib = terralib.includec("stdlib.h")
local stdio = terralib.includec("stdio.h")
local cstring = terralib.includec("string.h")
local c = regentlib.c

task model_write(step: int, sea_surface_now : region(ispace(int2d), uvt_field), sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
                 velocity_now : region(ispace(int2d), uv_field), grid: region(ispace(int2d), grid_fields)) 
     where reads(sea_surface_now.t, grid.{xt, yt}, sea_bed_to_mean_sea_level.t, velocity_now.{u,v}) do

--  var step : int32 = 4
  var filename : int8[50] 
  c.sprintf(&filename[0], "go2d_%05i_%05i.dat", step, 0)
  stdio.printf("%s\n",filename)
  var f = c.fopen(filename, 'w')
  var xlo = sea_surface_now.bounds.lo.x
  var xhi = sea_surface_now.bounds.hi.x
  var ylo = sea_surface_now.bounds.lo.y
  var yhi = sea_surface_now.bounds.hi.y

  for y=2, yhi+1 do
    for x=2, xhi+1 do
      var point : int2d = int2d({x,y})
      var rtmp1 = 0.5 * (velocity_now[point + {-1,0}].u + velocity_now[point].u) 
      var rtmp2 = 0.5 * (velocity_now[point + {0,-1}].v + velocity_now[point].v)
      c.fprintf(f, "%16.7e %16.7e %16.7e %16.7e %16.7e %16.7e\n", grid[point].xt, grid[point].yt, sea_bed_to_mean_sea_level[point].t, sea_surface_now[point].t, rtmp1, rtmp2)
    end
    cfprintf(f, "\n")
  end

  c.fclose(f)

end

--task main()
--regentlib.start(main)
