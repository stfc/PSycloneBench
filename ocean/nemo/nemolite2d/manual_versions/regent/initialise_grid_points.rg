import "regent"

require("model_init")

fspace uvt_field{
  u : double,
  v : double,
  t : double
}

fspace uvt_area_field{
  u : double,
  u_area: double,
  v : double,
  v_area : double,
  t : double,
  t_area : double
}

fspace uv_field{
  u : double,
  v : double,
}

fspace uv_area_field{
  u : double,
  u_area : double,
  v : double,
  v_area : double,
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


task setup_sea_surface_now( sea_surface_now: region(ispace(int2d), uvt_field), grid : region(ispace(int2d), grid_fields)) where reads(grid.{area_v, area_t, area_u}, sea_surface_now.t), writes(sea_surface_now.{u,v}) do 
  
  var full_partition = partition(equal, sea_surface_now, ispace(int2d,{1,1}))
  --Only set sea_surface values in the non-boundaries
  var centre_region = image(sea_surface_now, full_partition, calculate_internal_size)
  init_surface_now_u(centre_region[int2d({0,0})], grid)
  init_surface_now_v(centre_region[int2d({0,0})], grid)

end

task setup_sea_surface_after( sea_surface_after : region(ispace(int2d), uvt_field))

  -- Initialise values

end

task setup_sea_bed_to_mean_sea_level(  sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field))

  -- Initialise values

end

task setup_velocity_now( velocity_now : region(ispace(int2d), uv_field))

  -- Initialise values

end

task setup_velocity_after( velocity_after : region(ispace(int2d), uv_field))

  -- Initialise value

end
