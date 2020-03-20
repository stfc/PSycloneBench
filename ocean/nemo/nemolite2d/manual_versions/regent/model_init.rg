import "regent"

require("read_namelist")

local stdlib = terralib.includec("stdlib.h")
local stdio = terralib.includec("stdio.h")

fspace tmask_fields{
  tmask : int1d
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

task dump_tmask( tmask: region(ispace(int2d), tmask_fields)) where reads(tmask.tmask) do

  var f = stdio.fopen("tmask_0.dat", 'w')
  var xlo = tmask.bounds.lo.x
  var xhi = tmask.bounds.hi.x
  var ylo = tmask.bounds.lo.y
  var yhi = tmask.bounds.hi.y

  for x=xlo,xhi+1 do
    for y=ylo,yhi+1 do
      stdio.fprintf(f,"%i %i %i\n", x, y, tmask[int2d({x,y})].tmask)
    end
  end
  

end


task init_centre( tmask_centre : region(ispace(int2d), tmask_fields)) where writes( tmask_centre.tmask) do

--  __demand(__vectorize)
  for point in tmask_centre do
    tmask_centre[point].tmask = WET
  end
end

task init_centre_launcher( tmask_centre: region(ispace(int2d), tmask_fields)) where writes(tmask_centre.tmask) do

  var full_partition = partition(equal, tmask_centre, ispace(int2d, {4,4}))
  for point in ispace(int2d, {4,4}) do
    init_centre(full_partition[point])
  end

end

task init_west( tmask_west : region(ispace(int2d), tmask_fields)) where writes(tmask_west.tmask) do

  for point in tmask_west do
    tmask_west[point].tmask = SOLID_BOUNDARY
  end

end

task init_east( tmask_east : region(ispace(int2d), tmask_fields)) where writes(tmask_east.tmask) do

  for point in tmask_east do
    tmask_east[point].tmask = SOLID_BOUNDARY
  end
end

task init_north( tmask_north : region(ispace(int2d), tmask_fields)) where writes(tmask_north.tmask) do

  for point in tmask_north do
    tmask_north[point].tmask = SOLID_BOUNDARY
  end
end
task init_south( tmask_south : region(ispace(int2d), tmask_fields)) where writes(tmask_south.tmask) do

  for point in tmask_south do
    tmask_south[point].tmask = TIDAL_MASK
  end
end


task print_tmask( tmask_centre : region(ispace(int2d), tmask_fields) ) where reads( tmask_centre.tmask) do
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
  return rect2d({ private_bounds.lo + {1,1}, private_bounds.hi - {1,1} })
end
 
task calculate_west_boundary( private_bounds: rect2d) : rect2d
  return rect2d( { {private_bounds.lo.x, private_bounds.lo.y+1}, {private_bounds.lo.x, private_bounds.hi.y}})
end
task calculate_east_boundary( private_bounds: rect2d) : rect2d
  return rect2d( { {private_bounds.hi.x, private_bounds.lo.y+1}, {private_bounds.hi.x, private_bounds.hi.y}})
end
task calculate_south_boundary( private_bounds: rect2d) : rect2d
  return rect2d( { {private_bounds.lo.x, private_bounds.lo.y}, {private_bounds.hi.x, private_bounds.lo.y} })
end
task calculate_north_boundary( private_bounds: rect2d) : rect2d
  return rect2d( { {private_bounds.lo.x+1, private_bounds.hi.y}, {private_bounds.hi.x-1, private_bounds.hi.y} })
end

task main()

  var setup_data_is = ispace(int1d, 1)
  var setup_data = region(setup_data_is, setup_type)
  fill(setup_data.{jpiglo, jpjglo, nit000, nitend, record, jphgr_msh}, 0)
  fill(setup_data.{dx, dy, dep_const, rdt, cbfr, visc}, 0)
  init(setup_data)



--  var xdim : int32 = 5--setup_data[0].jpiglo
--  var ydim : int32 = 5--setup_data[0].jpjglo
  var xdim : int32 = setup_data[0].jpiglo
  var ydim : int32 = setup_data[0].jpjglo

  var grid_space = ispace(int2d, {x = xdim+2, y = ydim+2}, {x = 1, y = 1})
  var tmasks = region(grid_space, tmask_fields)
  var full_partition = partition(equal, tmasks, ispace(int2d, {1,1}))

  fill(tmasks.tmask, -2)

--  var partitioned_region = partition(equal, tmasks, ispace(int2d, {4,4}))
  var centre_region = image(tmasks, full_partition, calculate_internal_size)
  var west_region = image(tmasks, full_partition, calculate_west_boundary) 
  var east_region = image(tmasks, full_partition, calculate_east_boundary)
  var north_region = image(tmasks, full_partition, calculate_north_boundary)
  var south_region = image(tmasks, full_partition, calculate_south_boundary)
--  var centre_region_partitioned = partition(equal, tmasks, ispace(int2d,{4,4}))
--  var test = centre_region & centre_region_partitioned
--  for val in test.colors do
--  stdio.printf("%i %i\n", val.x, val.y)
--  end
--  for point in centre_region do
    init_centre_launcher(centre_region[int2d({0,0})])
    init_west(west_region[int2d({0,0})])
    init_east(east_region[int2d({0,0})])
    init_north(north_region[int2d({0,0})])
    init_south(south_region[int2d({0,0})])
--  end
--    print_tmask(tmasks)
    dump_tmask(tmasks)
 
end 

regentlib.start(main)
