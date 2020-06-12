import "regent"



require("bc_flather")
require("calculate_sea_surface_height_t")
require("calculate_velocity_fields")
require("initialise_grid_points")
require("model_init")
require("model_write")
require("read_config")
require("read_namelist")
require("update_sea_surface_t")
require("update_values")
require("update_velocity_boundary")

local c = regentlib.c
local pi = 3.1415926535897932
local g = 9.80665 
local omega = 7.292116e-05
local d2r = pi / 180.0


terra set_mappers()

end

task get_time() : double

  return c.legion_get_current_time_in_micros()
end

--INDEXING is ARRAY(X:M, Y:N)


task calculate_1_to_N_1_to_M(private_bounds : rect2d) : rect2d
  return rect2d({{1,1}, {private_bounds.hi.x-1, private_bounds.hi.y-1}})

end

task calculate_2_to_N_2_to_M( private_bounds : rect2d) : rect2d

  return rect2d({{2,2},  {private_bounds.hi.x-1, private_bounds.hi.y-1}})

end

task calculate_2_to_N_2_to_M1( private_bounds : rect2d) : rect2d

  return rect2d({{2,2},  {private_bounds.hi.x-2, private_bounds.hi.y-1}})

end

task calculate_2_to_N1_2_to_M( private_bounds : rect2d ) : rect2d


    return rect2d({{2,2}, {private_bounds.hi.x-1, private_bounds.hi.y-2}})
end

task calculate_1_to_full_1_to_M( private_bounds : rect2d) : rect2d

    return rect2d({{1,1}, {private_bounds.hi.x-1, private_bounds.hi.y}})

end

task calculate_1_to_N_1_to_full( private_bounds : rect2d) : rect2d

    return rect2d({{1,1}, {private_bounds.hi.x, private_bounds.hi.y-1}})

end

task calculate_halo_size( private_bounds: rect2d) : rect2d
  return rect2d({ private_bounds.lo - {1,1}, private_bounds.hi + {1,1} })
end

task flather_v_bounds( private_bounds : rect2d) : rect2d
  var suby = -1
  if private_bounds.lo.y == 1 then
    suby = 0
  end
  return rect2d({ private_bounds.lo - {0,suby}, private_bounds.hi + {0,1}})
end

task flather_u_bounds( private_bounds : rect2d) : rect2d
  var subx = -1
  if private_bounds.lo.x == 1 then
    subx = 0
  end
  return rect2d({ private_bounds.lo - {subx,0}, private_bounds.hi + {1,0}})

end

task main() 
  var conf : config = read_config()
  var single_point = ispace(int1d, 1)
  var setup_data = region(single_point, setup_type)
  fill(setup_data.{jpiglo, jpjglo, nit000, nitend, record, jphgr_msh}, 0)
  fill(setup_data.{dx, dy, dep_const, rdt, cbfr, visc}, 0)
  read_namelist(setup_data)

   var grid_space = ispace(int2d, {x = setup_data[0].jpiglo + 3,
                                   y = setup_data[0].jpjglo + 3},
                                   {x = 1, y = 1} )
   var grid = region(grid_space, grid_fields)
   var loop_conditions_data = region(grid_space, loop_conditions) 
--Initialise values (some of these are not changed after)
   fill(grid.tmask, -2)
   fill(grid.{dx_t, dx_u, dx_v}, setup_data[0].dx)
   fill(grid.{dy_t, dy_u, dy_v}, setup_data[0].dy)
   fill(grid.{area_t, area_u, area_v}, 0)
   fill(grid.{gphi_u, gphi_v}, 50.0)
   fill(grid.{xt, yt}, 0.0)
 --Initialise loop condition data - 1 means we don't use this value.
   fill(loop_conditions_data.{compute_vel_ufield,compute_vel_vfield}, int1d(1))
   fill(loop_conditions_data.{update_sea_surface_t}, int1d(1))
   fill(loop_conditions_data.{update_uvel_boundary, update_vvel_boundary}, int1d(1))
   fill(loop_conditions_data.{bc_flather_u, bc_flather_v}, int1d(1))
   fill(loop_conditions_data.{update_u_height, update_v_height}, int1d(1))


  --Initialise model
  model_init(grid, loop_conditions_data)

 var sea_surface = region(grid_space, uvt_time_field)
 fill(sea_surface.{u_now, u_after, v_now, v_after, t_now, t_after}, 0.0)
 setup_sea_surface(sea_surface, grid)

 var velocity = region(grid_space, uv_time_field)
 fill(velocity.{u_now, u_after, v_now, v_after}, 0.0)


  var sea_bed_to_mean_sea_level = region(grid_space, uvt_field)
  fill(sea_bed_to_mean_sea_level.{u, v, t}, setup_data[0].dep_const)
  setup_sea_bed_to_mean_sea_level( sea_bed_to_mean_sea_level )


--  c.printf("%i\n", __raw(velocity).tree_id)


  --Create the partitions we need for the various computations
  --Create the full partitions
  var full_ispace =  ispace(int2d, {x=1, y=1})
  var full_grid = partition(equal, grid, full_ispace)
  var full_sea_surface = partition(equal, sea_surface, full_ispace)
  var full_sea_bed_to_mean_sea_level = partition(equal, sea_bed_to_mean_sea_level, full_ispace)
  var full_velocity = partition(equal, velocity, full_ispace)

  --Create the 1 to N, 1 to M partitions
   var _1N1M_velocity = image(disjoint, incomplete, velocity, full_velocity, calculate_1_to_N_1_to_M)

  --Create the 2 to N 2 to M partitions
   var _2N2M_sea_surface = image(disjoint, incomplete, sea_surface, full_sea_surface, calculate_2_to_N_2_to_M)

  --Create the 2 to N 2 to M-1 partitions
   var _2N2M1_sea_surface = image(disjoint, incomplete, sea_surface, full_sea_surface, calculate_2_to_N_2_to_M1)
   var _2N2M1_velocity = image(disjoint, incomplete, velocity, full_velocity, calculate_2_to_N_2_to_M1)

  --Create the 2 to N-1 2 to M partitions
   var _2N12M_sea_surface = image(disjoint, incomplete, sea_surface, full_sea_surface, calculate_2_to_N1_2_to_M)
   var _2N12M_velocity = image(disjoint, incomplete, velocity, full_velocity, calculate_2_to_N1_2_to_M)

  --Create the 1 to N+1 1 to M partition
   var _FN1M_velocity = image(disjoint, incomplete, velocity, full_velocity, calculate_1_to_full_1_to_M)

  --Create the 1 to N 1 to M+1 partition
   var _1NFM_velocity = image(disjoint, incomplete, velocity, full_velocity, calculate_1_to_N_1_to_full)

  var tilesize : int
  tilesize = conf.t1
  c.printf("Tilesize 1: %i\n", tilesize)
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
  var _2N2M1_velocity_halos = image(disjoint, incomplete, velocity, partitioned_2N2M1_velocity, calculate_halo_size)
  var partitioned_2N12M_velocity = partition(equal, _2N12M_velocity[int2d({0,0})], partition_space)
  var _2N12M_velocity_halos = image(disjoint, incomplete, velocity, partitioned_2N12M_velocity, calculate_halo_size)

  tilesize = conf.t2
  c.printf("Tilesize 2: %i\n", tilesize)
  local_x = setup_data[0].jpiglo / tilesize
  if(local_x < 1) then
    local_x = 1
  end
  local_y = setup_data[0].jpjglo / tilesize
  if(local_y < 1) then
    local_y = 1
  end
  var partition_space2 = ispace(int2d, {x = local_x, y = local_y})
  var partitioned_2N2M_sea_surface = partition(equal, _2N2M_sea_surface[int2d({0,0})], partition_space2)
  var _2N2M_sea_surface_halos = image(disjoint, complete, sea_surface, partitioned_2N2M_sea_surface, calculate_halo_size)
  var partitioned_full_velocity = partition(equal, full_velocity[int2d({0,0})], partition_space2)
  var partitioned_full_sea_surface = partition(equal, full_sea_surface[int2d({0,0})], partition_space2)
  var partitioned_2N2M1_sea_surface = partition(equal, _2N2M1_sea_surface[int2d({0,0})], partition_space2)
  var _2N2M1_sea_surface_halos = image(disjoint, incomplete, sea_surface, partitioned_2N2M1_sea_surface, calculate_halo_size)
  var partitioned_2N12M_sea_surface = partition(equal, _2N12M_sea_surface[int2d({0,0})], partition_space2)
  var _2N12M_sea_surface_halos = image(disjoint, incomplete, sea_surface, partitioned_2N12M_sea_surface, calculate_halo_size)
 
  var partition_space3 = ispace(int2d, {x = 1, y = local_y})
  var partitioned_1NFM_velocity = partition(equal, _1NFM_velocity[int2d({0,0})], partition_space3)
  var flather_1NFM_velocity = image(disjoint, complete, velocity, partitioned_1NFM_velocity, flather_v_bounds)


  var partition_space4 = ispace(int2d, {x = local_x, y = 1})
  var partitioned_FN1M_velocity = partition(equal, _FN1M_velocity[int2d({0,0})], partition_space4)
  var flather_FN1M_velocity = image(disjoint, complete, velocity, partitioned_FN1M_velocity, flather_u_bounds)

  model_write( 0, sea_surface, sea_bed_to_mean_sea_level, velocity, grid, 0)
  var visc = setup_data[0].visc 
  var rdt  = setup_data[0].rdt
  var cbfr = setup_data[0].cbfr
  var point : int2d = int2d({0,0})

  -- update_velocity_ufield partition
  var loop_vel_ufield_partition = partition( loop_conditions_data.compute_vel_ufield, ispace(int1d, 1 ) )
  var vel_ufield_partitions = cross_product(partitioned_2N2M1_velocity, loop_vel_ufield_partition)
 
  var loop_vel_vfield_partition = partition( loop_conditions_data.compute_vel_vfield, ispace(int1d, 1) )
  var vel_vfield_partitions = cross_product(partitioned_2N12M_velocity, loop_vel_vfield_partition)

  --update sea_surface_t partition
  var loop_update_sst_partition = partition( loop_conditions_data.update_sea_surface_t, ispace(int1d, 1) )
  var update_sst_partitions = cross_product(partitioned_2N2M_sea_surface, loop_update_sst_partition )

  --update_vvel_boundary partition
  var loop_update_vvel_boundary_partition = partition( loop_conditions_data.update_vvel_boundary, ispace(int1d, 1) )
  var update_vvel_bound_partitions = cross_product( partitioned_1NFM_velocity, loop_update_vvel_boundary_partition)

  --update_uvel_boundary partition
  var loop_update_uvel_boundary_partition = partition( loop_conditions_data.update_uvel_boundary, ispace(int1d, 1) )
  var update_uvel_bound_partitions  = cross_product( partitioned_FN1M_velocity, loop_update_uvel_boundary_partition)

  --bc_flather_v partition
  var loop_update_bc_flather_v_partition = partition( loop_conditions_data.bc_flather_v, ispace(int1d, 1) )
  var bc_flather_v_partitions = cross_product( flather_1NFM_velocity , loop_update_bc_flather_v_partition)

  --bc_flather_u partition
  var loop_update_bc_flather_u_partition = partition( loop_conditions_data.bc_flather_u, ispace(int1d, 1) )
  var bc_flather_u_partitions = cross_product( flather_FN1M_velocity, loop_update_bc_flather_u_partition)

  --update_u_height partition
--  var loop_update_u_height_partition = partition( loop_conditions_data.update_u_height, ispace(int1d, 1) )
--  var update_u_height_partitions = cross_product( partitioned_2N2M1_sea_surface, loop_update_u_height_partition)

  --update_v_height partition
--  var loop_update_v_height_partition = partition( loop_conditions_data.update_v_height, ispace(int1d, 1) )
--  var update_v_height_partitions = cross_product( partitioned_2N12M_sea_surface, update_v_height_partitions )

  __fence(__execution, __block)
  var start_time = get_time()
  __fence(__execution, __block) 
  --Main timestepping loop to do!

  __demand(__trace) --, __spmd)
  for i = setup_data[0].nit000, setup_data[0].nitend+1 do
  
     __demand(__trace, __index_launch)
    for part in partition_space2 do
      calculate_sea_surface_t(partitioned_2N2M_sea_surface[part],
                             _2N2M_sea_surface_halos[part],
                              full_sea_bed_to_mean_sea_level[point],
                              full_velocity[point],
                              full_grid[point],
                              rdt)
                             
    end
--     __demand(__trace, __index_launch)
     for part in partition_space do
       update_velocity_ufield(vel_ufield_partitions[part][0],
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
--     __demand(__trace, __index_launch)
     for part in partition_space do 
       update_velocity_vfield(vel_vfield_partitions[part][0],
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
--    __fence(__execution, __block)
--    __demand(__trace , __index_launch)
    for part in partition_space2 do
    update_sea_surface_t(update_sst_partitions[part][0],
                         full_grid[point],
                         rdt,
                         i)
    end
    for part in partition_space3 do
    update_vvel_boundary(update_vvel_bound_partitions[part][0],
                         full_grid[point])
    end
    for part in partition_space4 do
    update_uvel_boundary(update_uvel_bound_partitions[part][0],
                         full_grid[point])
    end
    for part in partition_space3 do
    bc_flather_v_loop(bc_flather_v_partitions[part][0],
                      full_sea_bed_to_mean_sea_level[point],
                      full_sea_surface[point],
                      full_grid[point],
                      g)
    end
    for part in partition_space4 do
    bc_flather_u_loop(bc_flather_u_partitions[part][0],
                      full_sea_bed_to_mean_sea_level[point],
                      full_sea_surface[point],
                      full_grid[point],
                      g)
    end
--     __demand(__trace, __index_launch)
--     for part in partition_space2 do
--     update_velocity_and_t_height(partitioned_full_velocity[part],
--                                 partitioned_full_sea_surface[part])
--    end
    copy(velocity.{u_after, v_after},velocity.{u_now, v_now})
    copy(sea_surface.{t_after}, sea_surface.{t_now})

    __demand(__trace, __index_launch)
    for part in partition_space2 do
    update_u_height_launcher(partitioned_2N2M1_sea_surface[part],
                             full_grid[point],
                             _2N2M1_sea_surface_halos[part])
    end
    __demand(__trace, __index_launch)
    for part in partition_space2 do
    update_v_height_launcher(partitioned_2N12M_sea_surface[part],
                             full_grid[point],
                             _2N12M_sea_surface_halos[part])
    end
--    c.legion_runtime_issue_execution_fence(__runtime(), __context()) 
--    if( i % setup_data[0].record == 0) then
        model_write( i, sea_surface, sea_bed_to_mean_sea_level, velocity, grid, i % setup_data[0].record)
--    end

  end
  __fence(__execution, __block)
  var finish_time = get_time()
  __fence(__execution, __block)
  var time_taken = finish_time - start_time
  c.printf("Runtime is %f seconds\n", double(time_taken) / 1000000.0)
end


if os.getenv('SAVEOBJ') == '1' then
  local root_dir = arg[0]:match(".*/") or "./"
  local out_dir = (os.getenv('OBJNAME') and os.getenv('OBJNAME'):match('.*/')) or root_dir
  local link_flags = terralib.newlist({"-L" .. out_dir, "-lm", "-lgfortran", "-lgocean2d_io_mod"})

--  if os.getenv('STANDALONE') == '1' then
--    os.execute('cp ' .. os.getenv('LG_RT_DIR') .. '/../bindings/regent/' ..
--        regentlib.binding_library .. ' ' .. out_dir)
--  end

  local exe = os.getenv('OBJNAME') or "nemolite2D"
  regentlib.saveobj(main, exe, "executable", set_mappers, link_flags)
else
  regentlib.start(main)
end
