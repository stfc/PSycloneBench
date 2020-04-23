import "regent"


local stdlib = terralib.includec("stdlib.h")
local stdio = terralib.includec("stdio.h")
local cstring = terralib.includec("string.h")
local c = regentlib.c

task model_write(step: int, sea_surface : region(ispace(int2d), uvt_time_field), sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
                 velocity : region(ispace(int2d), uv_time_field), grid: region(ispace(int2d), grid_fields), l_write : int) 
     where reads(sea_surface.t_now, grid.{xt, yt}, sea_bed_to_mean_sea_level.t, velocity.{u_now,v_now}) do

  if(l_write == 0) then
    var filename : int8[50]
    var f : &c.FILE 
    c.sprintf(&filename[0], "go2d_%05i_%05i.dat", step, 0)
    f = c.fopen(filename, 'w')
    var xlo = sea_surface.bounds.lo.x
    var xhi = sea_surface.bounds.hi.x
    var ylo = sea_surface.bounds.lo.y
    var yhi = sea_surface.bounds.hi.y
  
    for y=2, yhi-1 do
      for x=2, xhi-1 do
        var point : int2d = int2d({x,y})
        var rtmp1 = 0.5 * (velocity[point + {-1,0}].u_now + velocity[point].u_now) 
        var rtmp2 = 0.5 * (velocity[point + {0,-1}].v_now + velocity[point].v_now)
        c.fprintf(f, "%16.6e %16.6e %16.6e %16.6e %16.6e %16.6e\n", grid[point].xt, grid[point].yt, sea_bed_to_mean_sea_level[point].t, sea_surface[point].t_now, rtmp1, rtmp2)
      end
        c.fprintf(f, "\n")
    end
  
    c.fclose(f)
  end
end

--task main()
--regentlib.start(main)
