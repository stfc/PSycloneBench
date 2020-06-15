import "regent"

require("initialise_grid_points")

local c = regentlib.c

local abs = c.fabs

--local function make_checksum_task(field_name)
--
----  local sum = regentlib.newsymbol(double, "sum")
----  local nil_sum = rexpr sum = 0.0 end
----  local loop_stuff = rexpr sum = sum + abs(input_region.[field_name]) end
--  local input_region = regentlib.newsymbol(region(ispace(int2d), uv_time_field), "input_region")
----  local task_name = regentlib.newsymbol(string, "checksum_".."field_name")
--
--  local task t( [input_region] ) where reads(input_region.[field_name]) do
--     var sum = 0.0
--     for x in [input_region] do
--     sum = sum + abs(input_region[x].[field_name])
--     end
--     return sum
--  end
--  return t
--end
--
--local u_after = "u_after"
--local v_after = "v_after"
--
--make_checksum_task(u_after)
--make_checksum_task(v_after)

task checksum_task( input_field : region(ispace(int2d), double)) : double
where reads(input_field) do

  var sum : double = 0.0
  for x in input_field do
    sum = sum + input_field[x]
  end
  return sum
end
