import "regent"

require("initialise_grid_points")
require("model_init")


local c = regentlib.c

sqrt = c.sqrt


--This is the SEVENTH loop
--Writes over N to M+1
--Reads from TODO
--Original code claims this is not parallelisable
task bc_flather_v_loop(velocity : region(ispace(int2d), uv_time_field),
                       sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
                       sea_surface : region(ispace(int2d), uvt_time_field), 
                       grid_region : region(ispace(int2d), grid_fields),
                       g : double)
     where -- sea_bed_to_mean_sea_level * sea_surface,
     writes(velocity.v_after), 
     reads(sea_bed_to_mean_sea_level.v, grid_region.tmask, sea_surface.v_now,
           velocity.v_after) do
--     DO jj = 1, N, 1
--       DO ji = 1, M+1, 1
--!          call bc_flather_v_code(ji,jj, &
--!                                 va%data, hv%data, sshn_v%data, &
--!                                 sshn_v%grid%tmask)
--          IF(sshn_t%grid%tmask(ji,jj) + sshn_t%grid%tmask(ji,jj+1) <= -1) cycle
--
--          IF(sshn_t%grid%tmask(ji,jj) < 0) THEN
--             jiv = jj + 1
--             va%data(ji,jj) = va%data(ji,jiv) + SQRT(g/hv%data(ji,jj)) * &
--                  (sshn_v%data(ji,jj) - sshn_v%data(ji,jiv))
--          ELSE IF(sshn_t%grid%tmask(ji,jj+1) < 0) THEN
--             jiv = jj - 1
--             va%data(ji,jj) = va%data(ji,jiv) + SQRT(g/hv%data(ji,jj)) * &
--                  (sshn_v%data(ji,jj) - sshn_v%data(ji,jiv))
--          END IF
--
--       END DO
--    END DO


     var xmin = velocity.bounds.lo.x
     var xmax = velocity.bounds.hi.x
     var ymin = velocity.bounds.lo.y
     var ymax = velocity.bounds.hi.y

     for jj = xmin, xmax+1 do
         for ji = ymin, ymax+1 do
             var point = int2d({ji,jj})
             if(grid_region[point].tmask + grid_region[point + {0,1}].tmask > int1d(-1)) then

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

end 



--This is the EIGTH loop
--Writes over N+1 to M
--Reads from TODO
task bc_flather_u_loop(velocity : region(ispace(int2d), uv_time_field),
                       sea_bed_to_mean_sea_level : region(ispace(int2d), uvt_field),
                       sea_surface : region(ispace(int2d), uvt_time_field),
                       grid_region : region(ispace(int2d), grid_fields),
                       g : double)
     where -- sea_bed_to_mean_sea_level * sea_surface, 
     writes(velocity.u_after),
     reads(sea_bed_to_mean_sea_level.u, grid_region.tmask, sea_surface.u_now,
           velocity.u_after) do


--DO jj = 1, N+1, 1
--       DO ji = 1, M, 1
--!          call bc_flather_u_code(ji,jj, &
--!                                 ua%data, hu%data, sshn_u%data, &
--!                                 sshn_u%grid%tmask)
--          ! Check whether this point lies within the domain
--          if(sshn_t%grid%tmask(ji,jj) + sshn_t%grid%tmask(ji+1,jj) <= -1) cycle
--
--          if(sshn_t%grid%tmask(ji,jj) < 0) then
--             ! Read from column to the right (East) of us
--             jiu = ji + 1
--             ua%data(ji,jj) = ua%data(jiu,jj) + sqrt(g/hu%data(ji,jj))* &
--                  (sshn_u%data(ji,jj) - sshn_u%data(jiu,jj))
--          else if(sshn_t%grid%tmask(ji+1,jj )< 0) then
--             ! Read from column to the left of us
--             jiu = ji - 1
--             ua%data(ji,jj) = ua%data(jiu,jj) + sqrt(g/hu%data(ji,jj)) * &
--                  (sshn_u%data(ji,jj) - sshn_u%data(jiu,jj))
--          end if
--       END DO
--    END DO


     var xmin = velocity.bounds.lo.x
     var xmax = velocity.bounds.hi.x
     var ymin = velocity.bounds.lo.y
     var ymax = velocity.bounds.hi.y

    for jj = xmin, xmax+1 do
        for ji = ymin, ymax+1 do
            var point = int2d({ji,jj})        
            if( grid_region[point].tmask + grid_region[point + {1,0}].tmask > int1d(-1)) then
            
                if(grid_region[point].tmask < int1d(0)) then
                    velocity[point].u_after = velocity[point + {1,0}].u_after
                                            + sqrt(g / sea_bed_to_mean_sea_level[point].u)
                                            * (sea_surface[point].u_now
                                               - sea_surface[point + {1,0}].u_now)
                elseif( grid_region[point + {1,0}].tmask < int1d(0)) then
                    velocity[point].u_after = velocity[point + {-1,0}].u_after
                                            + sqrt(g / sea_bed_to_mean_sea_level[point].u)
                                            * (sea_surface[point].u_now
                                               - sea_surface[point + {-1,0}].u_now)
                end
            end
         end
    end
end
