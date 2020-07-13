import "regent"

require("initialise_grid_points")

local c = regentlib.c

local abs = c.fabs

--It is possible to do these checksum tasks without field specific inputs, however this disabled RDIR which could lose performance

--Checksum task for u_after
task checksum_task_u_after( input_field : region(ispace(int2d), uv_time_field)) : double
where reads(input_field.u_after) do

  var sum : double = 0.0
  for x in input_field do
    sum = sum + abs(input_field[x].u_after)
  end
  return sum
end

--Checksum task for v_after
task checksum_task_v_after( input_field : region(ispace(int2d), uv_time_field)) : double
where reads(input_field.v_after) do

  var sum : double = 0.0
  for x in input_field do
    sum = sum + abs(input_field[x].v_after)
  end
  return sum
end
