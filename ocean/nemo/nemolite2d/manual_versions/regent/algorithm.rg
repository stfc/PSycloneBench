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

task calculate_halo_size( private_bounds: rect2d) : rect2d
  return rect2d({ private_bounds.lo - {1,1}, private_bounds.hi + {1,1} })
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

 var sea_surface = region(grid_space, uvt_time_field)
 fill(sea_surface.{u_now, u_after, v_now, v_after, t_now, t_after}, 0.0)
 setup_sea_surface(sea_surface, grid)

 var velocity = region(grid_space, uv_time_field)
 fill(velocity.{u_now, u_after, v_now, v_after}, 0.0)


  var sea_bed_to_mean_sea_level = region(grid_space, uvt_field)
  fill(sea_bed_to_mean_sea_level.{u, v, t}, setup_data[0].dep_const)
  setup_sea_bed_to_mean_sea_level( sea_bed_to_mean_sea_level )


--  c.printf("%i\n", __raw(velocity).tree_id)
  model_write( 0, sea_surface, sea_bed_to_mean_sea_level, velocity, grid, 0)


  --Create the partitions we need for the various computations
  --Create the full partitions
  var full_ispace =  ispace(int2d, {x=1, y=1})
  var full_grid = partition(equal, grid, full_ispace)
  var full_sea_surface = partition(equal, sea_surface, full_ispace)
  var full_sea_bed_to_mean_sea_level = partition(equal, sea_bed_to_mean_sea_level, full_ispace)
  var full_velocity = partition(equal, velocity, full_ispace)

  --Create the 1 to N, 1 to M partitions
   var _1N1M_velocity = image(velocity, full_velocity, calculate_1_to_N_1_to_M)

  --Create the 2 to N 2 to M partitions
   var _2N2M_sea_surface = image(sea_surface, full_sea_surface, calculate_2_to_N_2_to_M)

  --Create the 2 to N 2 to M-1 partitions
   var _2N2M1_sea_surface = image(sea_surface, full_sea_surface, calculate_2_to_N_2_to_M1)
   var _2N2M1_velocity = image(velocity, full_velocity, calculate_2_to_N_2_to_M1)

  --Create the 2 to N-1 2 to M partitions
   var _2N12M_sea_surface = image(sea_surface, full_sea_surface, calculate_2_to_N1_2_to_M)
   var _2N12M_velocity = image(velocity, full_velocity, calculate_2_to_N1_2_to_M)

  --Create the 1 to N+1 1 to M partition
   var _FN1M_velocity = image(velocity, full_velocity, calculate_1_to_full_1_to_M)

  --Create the 1 to N 1 to M+1 partition
   var _1NFM_velocity = image(velocity, full_velocity, calculate_1_to_N_1_to_full)

  var tilesize : int = 256 --128
  --Partition the velocity field as required for the update_velocity launcher
  var local_x : int = setup_data[0].jpiglo / tilesize
  if(local_x < 1) then
    local_x = 1
  end
  var local_y : int = setup_data[0].jpjglo / tilesize
  if(local_y < 1) then
    local_y = 1
  end
  var partition_space = ispace(int2d, {x = local_x, y = local_y})
  var partitioned_2N2M1_velocity = partition(equal, _2N2M1_velocity[int2d({0,0})], partition_space)
  var _2N2M1_velocity_halos = image(velocity, partitioned_2N2M1_velocity, calculate_halo_size)
  var partitioned_2N12M_velocity = partition(equal, _2N12M_velocity[int2d({0,0})], partition_space)
  var _2N12M_velocity_halos = image(velocity, partitioned_2N12M_velocity, calculate_halo_size)


  var visc : double = setup_data[0].visc 
  var rdt  : double = setup_data[0].rdt
  var cbfr : double = setup_data[0].cbfr
  var point : int2d = int2d({0,0})
  __fence(__execution, __block)
  var start_time = c.legion_get_current_time_in_micros()
  --Main timestepping loop to do!

  __demand(__trace)
  for i = setup_data[0].nit000, setup_data[0].nitend+1 do
  
     calculate_sea_surface_t(_2N2M_sea_surface[point],
                            full_sea_bed_to_mean_sea_level[point],
                            full_velocity[point],
                            full_grid[point],
                            setup_data[0].rdt,
                            setup_data[0].jpiglo,
                            setup_data[0].jpjglo )
     __demand(__trace, __index_launch)
     for part in partition_space do
       update_velocity_ufield(partitioned_2N2M1_velocity[part],
                              full_grid[point],
                              _2N2M1_velocity_halos[part],
                              full_sea_bed_to_mean_sea_level[point],
                              full_sea_surface[point],
                              visc,
                              omega,
                              g,
                              rdt,
                              cbfr,
                              d2r)
     end
     __demand(__trace, __index_launch)
     for part in partition_space do 
       update_velocity_vfield(partitioned_2N12M_velocity[part],
                              full_grid[point],
                              _2N12M_velocity_halos[part],
                              full_sea_bed_to_mean_sea_level[point],
                              full_sea_surface[point],
                              visc,
                              omega,
                              g,
                              rdt,
                              cbfr,
                              d2r)
                              
     end
--     update_velocity_ufield(_2N2M1_velocity[point],
--                           full_grid[point],
--                           full_velocity[point],
--                           full_sea_bed_to_mean_sea_level[point],
--                           full_sea_surface[point],
--                           visc,
--                           omega,
--                           g,
--                           rdt,
--                           cbfr,
--                           d2r)
--    update_velocity_vfield(_2N12M_velocity[point],
--                           full_grid[point],
--                           full_velocity[point],
--                           full_sea_bed_to_mean_sea_level[point],
--                           full_sea_surface[point],
--                           setup_data[0].visc,
--                           omega,
--                           g,
--                           setup_data[0].rdt,
--                           setup_data[0].cbfr,
--                           d2r)
    update_sea_surface_t(_2N2M_sea_surface[point],
                         full_grid[point],
                         setup_data[0].rdt,
                         i)
    update_uvel_boundary(_FN1M_velocity[point],
                         full_grid[point])
    update_vvel_boundary(_1NFM_velocity[point],
                         full_grid[point])
    bc_flather_v_loop(_1NFM_velocity[point],
                      full_sea_bed_to_mean_sea_level[point],
                      full_sea_surface[point],
                      full_grid[point],
                      g)
    bc_flather_u_loop(_FN1M_velocity[point],
                      full_sea_bed_to_mean_sea_level[point],
                      full_sea_surface[point],
                      full_grid[point],
                      g)
    update_velocity_and_t_height(full_velocity[point],
                                 full_sea_surface[point])
    update_u_height_launcher(_2N2M1_sea_surface[point],
                             full_grid[point])
    update_v_height_launcher(_2N12M_sea_surface[point],
                             full_grid[point])
--    c.legion_runtime_issue_execution_fence(__runtime(), __context()) 
--    if( i % setup_data[0].record == 0) then
        model_write( i, sea_surface, sea_bed_to_mean_sea_level, velocity, grid, i % setup_data[0].record)
--    end

  end
  __fence(__execution, __block)
  var finish_time = c.legion_get_current_time_in_micros()
  var time_taken = finish_time - start_time
  c.printf("Runtime is %f seconds\n", double(time_taken) / 1000000.0)
end



regentlib.start(main)
--regentlib.profile(main)
