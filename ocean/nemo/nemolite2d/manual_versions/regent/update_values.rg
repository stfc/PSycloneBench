import "regent"

require("initialise_grid_points")
require("model_init")

local c = regentlib.c

------
------ This requires the full field, 1 to N+1 and 1 to M+1
------ This is the NINTH loop
------
__demand(__leaf)
task update_velocity_and_t_height( velocity : region(ispace(int2d), uv_time_field), sea_surface : region(ispace(int2d), uvt_time_field)) where 
                                   writes(velocity.{u_now,v_now}, sea_surface.t_now),
                                   reads(velocity.{u_after,v_after}, sea_surface.t_after) do

--     do jj = 1, N+1, 1
--       do ji = 1, M+1, 1
--          un%data(ji,jj) = ua%data(ji,jj)
--          vn%data(ji,jj) = va%data(ji,jj)
--          sshn_t%data(ji,jj) = ssha%data(ji,jj)
--       end do
--    end do

  --TODO Create launcher function.
  __demand(__openmp, __vectorize)
  for point in velocity do
    velocity[point].u_now = velocity[point].u_after
    velocity[point].v_now = velocity[point].v_after
    sea_surface[point].t_now = sea_surface[point].t_after
  end
end

task get_u_height_launcher_bounds( private_bounds : rect2d) : rect2d
  return rect2d( {{2,2}, {private_bounds.hi.x, private_bounds.hi.y-1}})
end


--------
-------- This writes to a sub field, 2 to N, 2 to M-1
-------- This reads from the full field
--This is the TENTH loop
__demand(__leaf)
task update_u_height_launcher( sea_surface : region(ispace(int2d), uvt_time_field), grid_region : region(ispace(int2d), grid_fields) )
                     where writes(sea_surface.u_now), reads(grid_region.area_u, grid_region.area_t, grid_region.tmask, sea_surface.t_now) do
 --    do jj = 2, N, 1
--!dir$ vector always
--      do ji = 2, M-1, 1
--
--
--         if(sshn_t%grid%tmask(ji,jj) + &
--            sshn_t%grid%tmask(ji+1,jj) <= 0) cycle !jump over non-computational domain
--   var full_partition = partition(equal, sea_surface_now, ispace(int2d, {1,1}))
--   var u_height_region = image(sea_surface_now, full_partition, get_u_height_launcher_bounds)
--
--   var tmask_field_partition = partition(grid_region.tmask, ispace(int1d, {2, -1}))
__demand(__openmp)
  for point in sea_surface do
    if( grid_region[point].tmask + grid_region[point+{1,0}].tmask > int1d(0)) then
      --         IF(sshn_t%grid%tmask(ji,jj) * sshn_t%grid%tmask(ji+1,jj) > 0) THENi
--            rtmp1 = sshn_t%grid%area_t(ji,jj) * sshn_t%data(ji,jj) + &
--                 sshn_t%grid%area_t(ji+1,jj) * sshn_t%data(ji+1,jj)
--            sshn_u%data(ji,jj) = 0.5_go_wp * rtmp1 / sshn_t%grid%area_u(ji,jj)
--         ELSE IF(sshn_t%grid%tmask(ji,jj) <= 0) THEN
--            sshn_u%data(ji,jj) = sshn_t%data(ji+1,jj)
--         ELSE IF(sshn_t%grid%tmask(ji+1,jj) <= 0) THEN
--            sshn_u%data(ji,jj) = sshn_t%data(ji,jj)
--         END IF

        if(grid_region[point].tmask * grid_region[point + {1,0}].tmask > int1d(0)) then
           var rtmp1 : double = grid_region[point].area_t * sea_surface[point].t_now
                          + grid_region[point+{1,0}].area_t * sea_surface[point+{1,0}].t_now
           sea_surface[point].u_now = 0.5 * rtmp1 / grid_region[point].area_u
        elseif( grid_region[point].tmask <= int1d(0) ) then
           sea_surface[point].u_now = sea_surface[point + {1,0}].t_now
        elseif( grid_region[point + {1,0}].tmask <= int1d(0)) then
           sea_surface[point].u_now = sea_surface[point].t_now
        end


  
    end
  end


end


--------
-------- This writes to a sub field, 2 to N-1, 2 to M
-------- This reads from the full field
--This is the ELEVENTH loop
__demand(__leaf)
task update_v_height_launcher( sea_surface : region(ispace(int2d), uvt_time_field), grid_region : region(ispace(int2d), grid_fields) )
                     where writes(sea_surface.v_now), reads(grid_region.area_v, grid_region.area_t, grid_region.tmask, sea_surface.t_now) do


--    do jj = 2, N-1, 1
--!dir$ vector always
--      do ji = 2, M, 1
--
--!        call next_sshv_code(ji, jj,                   &
--!                            sshn_v%data, sshn_t%data, &
--!                            sshn_t%grid%tmask,        &
--!                            sshn_t%grid%area_t, sshn_t%grid%area_v)
--
--         if(sshn_t%grid%tmask(ji,jj) + &
--            sshn_t%grid%tmask(ji,jj+1) <= 0)  cycle !jump over non-computational domaini
            --TODO Hopefully this condition is handled through an "image" partition

__demand(__openmp)
  for point in sea_surface do
    if( grid_region[point].tmask + grid_region[point+{0,1}].tmask > int1d(0)) then
--         if(sshn_t%grid%tmask(ji,jj) * sshn_t%grid%tmask(ji,jj+1) > 0) then
--            rtmp1 = sshn_t%grid%area_t(ji,jj)*sshn_t%data(ji,jj) + &
--                 sshn_t%grid%area_t(ji,jj+1) * sshn_t%data(ji,jj+1)
--            sshn_v%data(ji,jj) = 0.5_go_wp * rtmp1 / sshn_t%grid%area_v(ji,jj)
--         else if(sshn_t%grid%tmask(ji,jj) <= 0) then
--            sshn_v%data(ji,jj) = sshn_t%data(ji,jj+1)
--         else if(sshn_t%grid%tmask(ji,jj+1) <= 0) then
--            sshn_v%data(ji,jj) = sshn_t%data(ji,jj)
--         end if
         if(grid_region[point].tmask * grid_region[point + {0,1}].tmask > int1d(0)) then
           var rtmp1 : double = grid_region[point].area_t * sea_surface[point].t_now
                          + grid_region[point+{0,1}].area_t * sea_surface[point+{0,1}].t_now
           sea_surface[point].v_now = 0.5 * rtmp1 / grid_region[point].area_v
        elseif( grid_region[point].tmask <= int1d(0) ) then
           sea_surface[point].v_now = sea_surface[point + {0,1}].t_now
        elseif( grid_region[point + {0,1}].tmask <= int1d(0)) then
           sea_surface[point].v_now = sea_surface[point].t_now
        end 
    end

  end



end 
