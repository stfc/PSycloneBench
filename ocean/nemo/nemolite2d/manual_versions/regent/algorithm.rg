import "regent"



require("bc_flather")
require("calculate_sea_surface_height_t")
require("calculate_velocity_fields")
require("initialise_grid_points")
require("model_init")
require("model_write")
require("read_namelist")
require("update_sea_surface_t")
require("update_values")
require("update_velocity_boundary")

local c = regentlib.c
local pi = 3.1415926535897932
local g = 9.80665 
local omega = 7.292116e-05
local d2r = pi / 180.0


task calculate_1_to_N_1_to_M(private_bounds : rect2d) : rect2d
  return rect2d({{1,1}, {private_bounds.hi.x-1, private_bounds.hi.y-1}})

end

task calculate_2_to_N_2_to_M( private_bounds : rect2d) : rect2d

  return rect2d({{2,2},  {private_bounds.hi.x-1, private_bounds.hi.y-1}})

end

task calculate_2_to_N_2_to_M1( private_bounds : rect2d) : rect2d

  return rect2d({{2,2},  {private_bounds.hi.x-1, private_bounds.hi.y-2}})

end

task calculate_2_to_N1_2_to_M( private_bounds : rect2d ) : rect2d


    return rect2d({{2,2}, {private_bounds.hi.x-2, private_bounds.hi.y-1}})
end

task calculate_1_to_full_1_to_M( private_bounds : rect2d) : rect2d

    return rect2d({{1,1}, {private_bounds.hi.x, private_bounds.hi.y-1}})

end

task calculate_1_to_N_1_to_full( private_bounds : rect2d) : rect2d

    return rect2d({{1,1}, {private_bounds.hi.x-1, private_bounds.hi.y}})

end

task main() 

  var single_point = ispace(int1d, 1)
  var setup_data = region(single_point, setup_type)
  fill(setup_data.{jpiglo, jpjglo, nit000, nitend, record, jphgr_msh}, 0)
  fill(setup_data.{dx, dy, dep_const, rdt, cbfr, visc}, 0)
  read_namelist(setup_data)

   var grid_space = ispace(int2d, {x = setup_data[0].jpiglo + 3,
                                   y = setup_data[0].jpjglo + 3},
                                   {x = 1, y = 1} )
   var grid = region(grid_space, grid_fields)
--Initialise values (some of these are not changed after)
   fill(grid.tmask, -2)
   fill(grid.{dx_t, dx_u, dx_v}, setup_data[0].dx)
   fill(grid.{dy_t, dy_u, dy_v}, setup_data[0].dy)
   fill(grid.{area_t, area_u, area_v}, 0)
   fill(grid.{gphi_u, gphi_v}, 50.0)
   fill(grid.{xt, yt}, 0.0)

  --Initialise model
  model_init(grid)

  var sea_surface_now = region(grid_space, uvt_field)
  fill(sea_surface_now.{u, v, t}, 0.0) 
  setup_sea_surface_now( sea_surface_now, grid )

  var sea_surface_after = region(grid_space, uvt_field)
  fill(sea_surface_after.{u,v,t}, 0.0)
  setup_sea_surface_after( sea_surface_after )

  var sea_bed_to_mean_sea_level = region(grid_space, uvt_field)
  fill(sea_bed_to_mean_sea_level.{u, v, t}, setup_data[0].dep_const)
  setup_sea_bed_to_mean_sea_level( sea_bed_to_mean_sea_level )

  var velocity_now = region(grid_space, uv_field)
  fill(velocity_now.{u, v}, 0.0)
  setup_velocity_now( velocity_now )

  var velocity_after = region(grid_space, uv_field) 
  fill(velocity_after.{u, v}, 0.0)
  setup_velocity_after( velocity_after )

--  c.printf("%i\n", __raw(velocity_now).tree_id)
--  c.printf("%i\n", __raw(velocity_after).tree_id)
  model_write( 0, sea_surface_now, sea_bed_to_mean_sea_level, velocity_now, grid, 0)


  --Create the partitions we need for the various computations
  --Create the full partitions
  var full_ispace =  ispace(int2d, {x=1, y=1})
  var full_grid = partition(equal, grid, full_ispace)
  var full_sea_surface_now = partition(equal, sea_surface_now, full_ispace)
  var full_sea_surface_after = partition(equal, sea_surface_after, full_ispace)
  var full_sea_bed_to_mean_sea_level = partition(equal, sea_bed_to_mean_sea_level, full_ispace)
  var full_velocity_now = partition(equal, velocity_now, full_ispace)
  var full_velocity_after = partition(equal, velocity_after, full_ispace)

  --Create the 1 to N, 1 to M partitions
   var _1N1M_grid = image(grid, full_grid, calculate_1_to_N_1_to_M)
   var _1N1M_sea_surface_now = image(sea_surface_now, full_sea_surface_now, calculate_1_to_N_1_to_M)
   var _1N1M_sea_surface_after = image(sea_surface_after, full_sea_surface_after, calculate_1_to_N_1_to_M)
   var _1N1M_sea_bed_to_mean_sea_level = image(sea_bed_to_mean_sea_level, full_sea_bed_to_mean_sea_level, calculate_1_to_N_1_to_M)
   var _1N1M_velocity_now = image(velocity_now, full_velocity_now, calculate_1_to_N_1_to_M)
   var _1N1M_velocity_after = image(velocity_after, full_velocity_after, calculate_1_to_N_1_to_M)

  --Create the 2 to N 2 to M partitions
   var _2N2M_grid = image(grid, full_grid, calculate_2_to_N_2_to_M)
   var _2N2M_sea_surface_now = image(sea_surface_now, full_sea_surface_now, calculate_2_to_N_2_to_M)
   var _2N2M_sea_surface_after = image(sea_surface_after, full_sea_surface_after, calculate_2_to_N_2_to_M)
   var _2N2M_sea_bed_to_mean_sea_level = image(sea_bed_to_mean_sea_level, full_sea_bed_to_mean_sea_level, calculate_2_to_N_2_to_M)
   var _2N2M_velocity_now = image(velocity_now, full_velocity_now, calculate_2_to_N_2_to_M)
   var _2N2M_velocity_after = image(velocity_after, full_velocity_after, calculate_2_to_N_2_to_M)

  --Create the 2 to N 2 to M-1 partitions
   var _2N2M1_grid = image(grid, full_grid, calculate_2_to_N_2_to_M)
   var _2N2M1_sea_surface_now = image(sea_surface_now, full_sea_surface_now, calculate_2_to_N_2_to_M1)
   var _2N2M1_sea_surface_after = image(sea_surface_after, full_sea_surface_after, calculate_2_to_N_2_to_M1)
   var _2N2M1_sea_bed_to_mean_sea_level = image(sea_bed_to_mean_sea_level, full_sea_bed_to_mean_sea_level, calculate_2_to_N_2_to_M1)
   var _2N2M1_velocity_now = image(velocity_now, full_velocity_now, calculate_2_to_N_2_to_M1)
   var _2N2M1_velocity_after = image(velocity_after, full_velocity_after, calculate_2_to_N_2_to_M1)

  --Create the 2 to N-1 2 to M partitions
   var _2N12M_grid = image(grid, full_grid, calculate_2_to_N1_2_to_M)
   var _2N12M_sea_surface_now = image(sea_surface_now, full_sea_surface_now, calculate_2_to_N1_2_to_M)
   var _2N12M_sea_surface_after = image(sea_surface_after, full_sea_surface_after, calculate_2_to_N1_2_to_M)
   var _2N12M_sea_bed_to_mean_sea_level = image(sea_bed_to_mean_sea_level, full_sea_bed_to_mean_sea_level, calculate_2_to_N1_2_to_M)
   var _2N12M_velocity_now = image(velocity_now, full_velocity_now, calculate_2_to_N1_2_to_M)
   var _2N12M_velocity_after = image(velocity_after, full_velocity_after, calculate_2_to_N1_2_to_M)

  --Create the 1 to N+1 1 to M partition
   var _FN1M_grid = image(grid, full_grid, calculate_1_to_full_1_to_M)
   var _FN1M_sea_surface_now = image(sea_surface_now, full_sea_surface_now, calculate_1_to_full_1_to_M)
   var _FN1M_sea_surface_after = image(sea_surface_after, full_sea_surface_after, calculate_1_to_full_1_to_M)
   var _FN1M_sea_bed_to_mean_sea_level = image(sea_bed_to_mean_sea_level, full_sea_bed_to_mean_sea_level, calculate_1_to_full_1_to_M)
   var _FN1M_velocity_now = image(velocity_now, full_velocity_now, calculate_1_to_full_1_to_M)
   var _FN1M_velocity_after = image(velocity_after, full_velocity_after, calculate_1_to_full_1_to_M)

  --Create the 1 to N 1 to M+1 partition
   var _1NFM_grid = image(grid, full_grid, calculate_1_to_N_1_to_full)
   var _1NFM_sea_surface_now = image(sea_surface_now, full_sea_surface_now, calculate_1_to_N_1_to_full)
   var _1NFM_sea_surface_after = image(sea_surface_after, full_sea_surface_after, calculate_1_to_N_1_to_full)
   var _1NFM_sea_bed_to_mean_sea_level = image(sea_bed_to_mean_sea_level, full_sea_bed_to_mean_sea_level, calculate_1_to_N_1_to_full)
   var _1NFM_velocity_now = image(velocity_now, full_velocity_now, calculate_1_to_N_1_to_full)
   var _1NFM_velocity_after = image(velocity_after, full_velocity_after, calculate_1_to_N_1_to_full)

  
  var point : int2d = int2d({0,0})
  __fence(__execution, __block)
  var start_time = c.legion_get_current_time_in_micros()
  --Main timestepping loop to do!

  for i = setup_data[0].nit000, setup_data[0].nitend+1 do
  
     calculate_sea_surface_t(_2N2M_sea_surface_after[point],
                          _1N1M_sea_surface_now[point],
                          _1N1M_sea_bed_to_mean_sea_level[point],
                          _1N1M_velocity_now[point],
                          _1N1M_grid[point],
                          setup_data[0].rdt)
     update_velocity_ufield_launcher(_2N2M1_velocity_after[point],
                           full_grid[point],
                           full_velocity_now[point],
                           full_sea_bed_to_mean_sea_level[point],
                           full_sea_surface_now[point],
                           full_sea_surface_after[point],
                           setup_data[0].visc,
                           omega,
                           g,
                           setup_data[0].rdt,
                           setup_data[0].cbfr,
                           d2r)
    update_velocity_vfield_launcher(_2N12M_velocity_after[point],
                           full_grid[point],
                           full_velocity_now[point],
                           full_sea_bed_to_mean_sea_level[point],
                           full_sea_surface_now[point],
                           full_sea_surface_after[point],
                           setup_data[0].visc,
                           omega,
                           g,
                           setup_data[0].rdt,
                           setup_data[0].cbfr,
                           d2r)
    update_sea_surface_t(_2N2M_sea_surface_after[point],
                         full_grid[point],
                         setup_data[0].rdt,
                         i)
    update_uvel_boundary(_FN1M_velocity_after[point],
                         full_grid[point])
    update_vvel_boundary(_FN1M_velocity_after[point],
                         full_grid[point])
    bc_flather_v_loop(_1NFM_velocity_after[point],
                      full_sea_bed_to_mean_sea_level[point],
                      full_sea_surface_now[point],
                      full_grid[point],
                      g)
    bc_flather_u_loop(_FN1M_velocity_after[point],
                      full_sea_bed_to_mean_sea_level[point],
                      full_sea_surface_now[point],
                      full_grid[point],
                      g)
    update_velocity_and_t_height(full_velocity_now[point],
                                 full_velocity_after[point],
                                 full_sea_surface_after[point],
                                 full_sea_surface_now[point])
    update_u_height_launcher(_2N2M1_sea_surface_now[point],
                             full_grid[point])
    update_v_height_launcher(_2N12M_sea_surface_now[point],
                             full_grid[point])
--    c.legion_runtime_issue_execution_fence(__runtime(), __context()) 
    if( i % setup_data[0].record == 0) then
        model_write( i, sea_surface_now, sea_bed_to_mean_sea_level, velocity_now, grid, i % setup_data[0].record)
    end

  end
  __fence(__execution, __block)
  var finish_time = c.legion_get_current_time_in_micros()
  var time_taken = finish_time - start_time
  c.printf("Runtime is %f seconds\n", double(time_taken) / 1000000.0)
end


regentlib.start(main)
--regentlib.profile(main)
