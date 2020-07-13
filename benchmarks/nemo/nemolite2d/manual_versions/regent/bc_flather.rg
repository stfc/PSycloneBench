import "regent"

require("initialise_grid_points")
require("model_init")



sqrt = regentlib.sqrt(double)


--This is the SEVENTH loop
__demand(__leaf)
task bc_flather_v_loop(velocity : region(ispace(int2d), uv_time_field),
                       sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
                       sea_surface : region(ispace(int2d), uvt_time_field), 
                       grid_region : region(ispace(int2d), grid_fields),
                       g : double)
     where
     writes(velocity.v_after), 
     reads(sea_bed_to_mean_sea_level.v, grid_region.tmask, sea_surface.v_now,
           velocity.v_after) do

--   Find the bounds of this partition to loop over
     var xmin = velocity.bounds.lo.x
     var xmax = velocity.bounds.hi.x
     var ymin = velocity.bounds.lo.y
     var ymax = velocity.bounds.hi.y-1
     for jj = ymin, ymax do
         for ji = xmin, xmax do
             var point = int2d({ji,jj})

                if(grid_region[point].tmask < int1d(0)) then
                  velocity[point].v_after = velocity[point + {0,1}].v_after
                                           + sqrt(g / sea_bed_to_mean_sea_level[point].v)
                                           * (sea_surface[point].v_now
                                              - sea_surface[point + {0,1}].v_now)
                elseif(grid_region[point + {0,1}].tmask < int1d(0) )then
                  velocity[point].v_after = velocity[point + {0,-1}].v_after
                                          + sqrt(g / sea_bed_to_mean_sea_level[point].v)
                                          * (sea_surface[point].v_now
                                            - sea_surface[point + {0,-1}].v_now)
                end

         end
     end

end 



--This is the EIGTH loop
__demand(__leaf)
task bc_flather_u_loop(velocity : region(ispace(int2d), uv_time_field),
                       sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
                       sea_surface : region(ispace(int2d), uvt_time_field),
                       grid_region : region(ispace(int2d), grid_fields),
                       g : double)
     where
     writes(velocity.u_after),
     reads(sea_bed_to_mean_sea_level.u, grid_region.tmask, sea_surface.u_now,
           velocity.u_after) do


--   Find the bounds of this partition to loop over
     var xmin = velocity.bounds.lo.x
     var xmax = velocity.bounds.hi.x-1
     var ymin = velocity.bounds.lo.y
     var ymax = velocity.bounds.hi.y
    for jj = ymin, ymax do
        for ji = xmin, xmax do
            var point = int2d({ji,jj})        

                if(grid_region[point].tmask < int1d(0)) then
--              Read from column to the right (East) of us
                    velocity[point].u_after = velocity[point + {1,0}].u_after
                                            + sqrt(g / sea_bed_to_mean_sea_level[point].u)
                                            * (sea_surface[point].u_now
                                               - sea_surface[point + {1,0}].u_now)
                elseif( grid_region[point + {1,0}].tmask < int1d(0)) then
--              Read from column to the left of us
                    velocity[point].u_after = velocity[point + {-1,0}].u_after
                                            + sqrt(g / sea_bed_to_mean_sea_level[point].u)
                                            * (sea_surface[point].u_now
                                               - sea_surface[point + {-1,0}].u_now)
                end
         end
    end
end
