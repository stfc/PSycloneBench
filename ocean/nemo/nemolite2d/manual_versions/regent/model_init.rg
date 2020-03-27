import "regent"

require("read_namelist")

local c = regentlib.c
local stdlib = terralib.includec("stdlib.h")
local stdio = terralib.includec("stdio.h")

fspace grid_fields{
  tmask : int1d,
  dx_t : double,
  dy_t : double,
  dx_u : double,
  dy_u : double,
  dx_v : double,
  dy_v : double,
  area_t : double,
  area_u : double, 
  area_v : double,
  gphi_u : double,
  gphi_v : double,
  xt : double,
  yt : double
-- gphi_f  : double  
--  dx_f : double -- Don't use f field in this code
--  dy_f : double
}

local SOLID_BOUNDARY = 0
local WET  = 1
local TIDAL_MASK = -1

--For this code x is East/West y is North/South
-- (1, N) --- --- -- -- (M,N)
--   |                    |
--   |                    |
--   |                    |
--   |                    |
--   |                    |
-- (1, 1) --- --- -- -- (M, 1)

task dump_tmask( tmask: region(ispace(int2d), grid_fields)) where reads(tmask.tmask) do

  var f = stdio.fopen("tmask_0.dat", 'w')
  var xlo = tmask.bounds.lo.x
  var xhi = tmask.bounds.hi.x
  var ylo = tmask.bounds.lo.y
  var yhi = tmask.bounds.hi.y

  for x=xlo,xhi do
    for y=ylo,yhi do
      stdio.fprintf(f,"%i %i %i\n", x, y, tmask[int2d({x,y})].tmask)
    end
  end
  

end

task init_grid_coordinates( tmask_region : region(ispace(int2d), grid_fields) ) where writes(tmask_region.{xt, yt}), reads(tmask_region.{dx_t, dy_t}) do

  for point in tmask_region do
     var xloc = (point.x-1) * tmask_region[point].dx_t
     var yloc = (point.y-1) * tmask_region[point].dy_t
     tmask_region[point].xt = xloc
     tmask_region[point].yt = yloc
  end

end

task init_grid_coordinates_launcher( tmask_full : region(ispace(int2d), grid_fields)) where writes(tmask_full.{xt, yt}), reads(tmask_full.{dx_t, dy_t}) do

  
  var full_partition = partition(equal, tmask_full, ispace(int2d, {4,4}))
  for point in ispace(int2d, {4,4}) do
    init_grid_coordinates(full_partition[point])
  end
end


task init_grid_areas( tmask_region : region(ispace(int2d), grid_fields) )  where writes(tmask_region.{area_t, area_u, area_v}),
   reads(tmask_region.{dx_t, dy_t, dx_u, dy_u, dx_v, dy_v}) do

  __demand(__vectorize)
  for point in tmask_region do
    tmask_region[point].area_t = tmask_region[point].dx_t * tmask_region[point].dy_t
    tmask_region[point].area_u = tmask_region[point].dx_u * tmask_region[point].dy_u
    tmask_region[point].area_v = tmask_region[point].dx_v * tmask_region[point].dy_v
  end 


end

task init_grid_areas_launcher( tmask_full : region(ispace(int2d), grid_fields)) where writes(tmask_full.{area_t, area_u, area_v}),
   reads(tmask_full.{dx_t, dy_t, dx_u, dy_u, dx_v, dy_v}) do

  var full_partition = partition(equal, tmask_full, ispace(int2d, {4,4}))
  for point in ispace(int2d, {4,4}) do
    init_grid_areas(full_partition[point])
  end
end

task init_centre( tmask_centre : region(ispace(int2d), grid_fields)) where writes( tmask_centre.tmask) do

--  __demand(__vectorize)
  for point in tmask_centre do
    tmask_centre[point].tmask = WET
  end
end

task init_centre_launcher( tmask_centre: region(ispace(int2d), grid_fields)) where writes(tmask_centre.tmask) do

  var full_partition = partition(equal, tmask_centre, ispace(int2d, {4,4}))
  for point in ispace(int2d, {4,4}) do
    init_centre(full_partition[point])
  end

end


task init_west( tmask_west : region(ispace(int2d), grid_fields)) where writes(tmask_west.tmask) do

  for point in tmask_west do
    tmask_west[point].tmask = SOLID_BOUNDARY
  end

end

task init_east( tmask_east : region(ispace(int2d), grid_fields)) where writes(tmask_east.tmask) do

  for point in tmask_east do
    tmask_east[point].tmask = SOLID_BOUNDARY
  end
end

task init_north( tmask_north : region(ispace(int2d), grid_fields)) where writes(tmask_north.tmask) do

  for point in tmask_north do
    tmask_north[point].tmask = SOLID_BOUNDARY
  end
end
task init_south( tmask_south : region(ispace(int2d), grid_fields)) where writes(tmask_south.tmask) do

  for point in tmask_south do
    tmask_south[point].tmask = TIDAL_MASK
  end
end


task print_tmask( tmask_centre : region(ispace(int2d), grid_fields) ) where reads( tmask_centre.tmask) do
  var xlo = tmask_centre.bounds.lo.x
  var xhi = tmask_centre.bounds.hi.x
  var ylo = tmask_centre.bounds.lo.y
  var yhi = tmask_centre.bounds.hi.y

  for y=yhi, ylo-1, -1 do 
    for x = xlo, xhi+1 do
      stdio.printf("%i ", tmask_centre[int2d({x,y})].tmask)
    end
    stdio.printf("\n")
  end
end

task calculate_internal_size( private_bounds: rect2d) : rect2d
  return rect2d({ private_bounds.lo + {1,1}, private_bounds.hi - {2,2} })
end
 
task calculate_west_boundary( private_bounds: rect2d) : rect2d
  return rect2d( { {private_bounds.lo.x, private_bounds.lo.y+1}, {private_bounds.lo.x, private_bounds.hi.y-1}})
end
task calculate_east_boundary( private_bounds: rect2d) : rect2d
  return rect2d( { {private_bounds.hi.x-1, private_bounds.lo.y+1}, {private_bounds.hi.x, private_bounds.hi.y-1}})
end
task calculate_south_boundary( private_bounds: rect2d) : rect2d
  return rect2d( { {private_bounds.lo.x, private_bounds.lo.y}, {private_bounds.hi.x, private_bounds.lo.y} })
end
task calculate_north_boundary( private_bounds: rect2d) : rect2d
  return rect2d( { {private_bounds.lo.x, private_bounds.hi.y-1}, {private_bounds.hi.x, private_bounds.hi.y} })
end

task model_init( grid : region(ispace(int2d), grid_fields)) where
    writes(grid.{tmask, dx_t, dx_u, dx_v, dy_t, dy_t, dy_v, gphi_u, gphi_v, xt, yt, area_t, area_u, area_v}), reads(grid.{tmask, dx_t, dx_u, dx_v, dy_t, dy_t, dy_u, dy_v, gphi_u, gphi_v, xt, yt, area_t ,area_u, area_v}) do

    var full_partition = partition(equal, grid, ispace(int2d, {1,1}))
  
    var centre_region = image(grid, full_partition, calculate_internal_size)
    var west_region = image(grid, full_partition, calculate_west_boundary)
    var east_region = image(grid, full_partition, calculate_east_boundary)
    var north_region = image(grid, full_partition, calculate_north_boundary)
    var south_region = image(grid, full_partition, calculate_south_boundary)
    init_centre_launcher(centre_region[int2d({0,0})])

    init_west(west_region[int2d({0,0})]) 

    init_east(east_region[int2d({0,0})])

    init_north(north_region[int2d({0,0})])

    init_south(south_region[int2d({0,0})])

    init_grid_coordinates(grid)
--    for point in ispace(int2d,{1,1}) do
--      init_grid_areas(full_partition[point])
--    end
--    init_grid_areas_launcher(full_partition[int2d({0,0})])
  __demand(__vectorize)
  for point in grid do
    grid[point].area_t = grid[point].dx_t * grid[point].dy_t
    grid[point].area_u = grid[point].dx_u * grid[point].dy_u
    grid[point].area_v = grid[point].dx_v * grid[point].dy_v
  end
  dump_tmask(grid)

end
