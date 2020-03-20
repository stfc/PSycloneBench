import "regent"

require("model_init")

fspace uvt_field{
  u : float
  v : float
  t : float
}

fspace uvt_area_field{
  u : float
  u_area: float
  v : float
  v_area : float
  t : float
  t_area : float
}

fspace uv_field{
  u : float
  v : float
}

fspace uv_area_field{
  u : float
  u_area : float
  v : float
  v_area : float
}

-- Import max/min for Terra
max = regentlib.fmax
min = regentlib.fmin


task calculate_internal_size( private_bounds: rect2d) : rect2d
  return rect2d({ private_bounds.lo + {1,1}, private_bounds.hi - {1,1} })
end




task init_surface_now_u( sea_surface_now : region(ispace(int2d), uvt_field), grid_region : region(ispace(int2d), grid_fields) ) where writes (sea_surface_now.u),
                                                                                                    reads(grid_region.{area_u, area_t}, sea_surface_now.t) do

--TODO partition this for parallelism
  for point in sea_surface_now do
    var rtmp = grid_region[point + {1,0}].area_t * sea_surface_now[point+{1,0}].t + grid_region[point].area_t * sea_surface_now[point].t
    sea_surface_now[point].u = 0.5 * rtmp / grid_region[point].area_u 
  end


end

task init_surface_now_v( sea_surface_now : region(ispace(int2d), uvt_field), grid_region : region(ispace(int2d), grid_fields) ) where writes (sea_surface_now.v),
                                                                                                    reads(grid_region.{area_v, area_t}, sea_surface_now.t) do

--TODO partition this for parallelism
  for point in sea_surface_now do
    var rtmp = grid_region[point + {0,1}].area_t * sea_surface_now[point+{0,1}].t + grid_region[point].area_t * sea_surface_now[point].t
    sea_surface_now[point].v = 0.5 * rtmp / grid_region[point].area_v 
  end


end


task setup_sea_surface_now( grid_space: ispace(int2d), region(ispace(int2d, grid_fields))) : region(ispace(int2d), uvt_field) 

  var sea_surface_now = region(grid_space, uvt_field)

  -- Initialise values
  fill(sea_surface_now.{u, v, t}, 0.0)
  --Only set sea_surface values in the non-boundaries
  var centre_region = image(tmasks, full_partition, calculate_internal_size)
  init_surface_now_u(centre_region[int2d({0,0})])
  return sea_surface_now

end

task setup_sea_surface_after( grid_space: ispace(int2d)) : region(ispace(int2d), uvt_field)

  var sea_surface_after = region(grid_space, uvt_field)
  -- Initialise values
  fill(sea_surface_after.{u,v,t}, 0.0)

  return sea_surface_after

end

task setup_sea_bed_to_mean_sea_level(  grid_space: ispace(int2d)) : region(ispace(int2d), uvt_field)

  var sea_bed_to_mean_sea_level =  region(grid_space, uvt_field)
  -- Initialise values
  fill(sea_bed_to_mean_sea_level.{u, v, t}, 0.0)

  return sea_bed_to_mean_sea_level
end

task setup_velocity_now( grid_space: ispace(int2d)) : region(ispace(int2d), uv_field)

  var velocity_now =  region(grid_space, uv_field)
  -- Initialise values
  fill(velocity_now.{u, v}, 0.0)

  return velocity_now
end

task setup_velocity_after( grid_space: ispace(int2d)) : region(ispace(int2d), uv_field)

  var velocity_after = region(grid_space, uv_field)
  -- Initialise value
  fill(velocity_after.{u, v}, 0.0)

  return velocity_after

end
