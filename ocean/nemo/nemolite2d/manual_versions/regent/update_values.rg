import "regent"

require("initialise_grid_points")
require("model_init")


------
------ This requires the full field, 1 to N+1 and 1 to M+1
------ This is the NINTH loop
------
task update_velocity_and_t_height( velocity_now : region(ispace(int2d), uv_field), velocity_after : region(ispace(int2d), uv_field), sea_surface_after : region(ispace(int2d), uvt_field),
                                   sea_surface_now :  region(ispace(int2d), uvt_field)) where writes(velocity_now.{u,v}, sea_surface_now.t),
                                   reads(velocity_after.{u,v}, sea_surface_after.t) do

--     do jj = 1, N+1, 1
--       do ji = 1, M+1, 1
--          un%data(ji,jj) = ua%data(ji,jj)
--          vn%data(ji,jj) = va%data(ji,jj)
--          sshn_t%data(ji,jj) = ssha%data(ji,jj)
--       end do
--    end do

  --TODO Create launcher function.
  for point in velocity_now do
    velocity_now[point].u = velocity_after[point].u
    velocity_now[point].v = velocity_after[point].v
    sea_surface_now[point].t = sea_surface_after[point].t
  end
end

task get_u_height_launcher_bounds( private_bounds : rect2d) : rect2d
  return rect2d( {{2,2}, {private_bounds.hi.x, private_bounds.hi.y-1}})
end


--------
-------- This writes to a sub field, 2 to N, 2 to M-1
-------- This reads from the full field
--This is the TENTH loop
task update_u_height_launcher( sea_surface_now : region(ispace(int2d), uvt_field), grid_region : region(ispace(int2d), grid_fields) )
                     where writes(sea_surface_now.u), reads(grid_region.area_u, grid_region.area_t, grid_region.tmask, sea_surface_now.t) do
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
  for point in sea_surface_now do
    if( grid_region[point].tmask + grid_region[point+{1,0}].tmask > 0) then
      --         IF(sshn_t%grid%tmask(ji,jj) * sshn_t%grid%tmask(ji+1,jj) > 0) THENi
--            rtmp1 = sshn_t%grid%area_t(ji,jj) * sshn_t%data(ji,jj) + &
--                 sshn_t%grid%area_t(ji+1,jj) * sshn_t%data(ji+1,jj)
--            sshn_u%data(ji,jj) = 0.5_go_wp * rtmp1 / sshn_t%grid%area_u(ji,jj)
--         ELSE IF(sshn_t%grid%tmask(ji,jj) <= 0) THEN
--            sshn_u%data(ji,jj) = sshn_t%data(ji+1,jj)
--         ELSE IF(sshn_t%grid%tmask(ji+1,jj) <= 0) THEN
--            sshn_u%data(ji,jj) = sshn_t%data(ji,jj)
--         END IF

        if(grid_region[point].tmask * grid_region[point + {1,0}].tmask > 0) then
           var rtmp1 : float = grid_region[point].area_t * sea_surface_now[point].t
                          + grid_region[point+{1,0}].area_t * sea_surface_now[point+{1,0}.t
           sea_surface_now[point].u = 0.5 * rtmp1 / grid_region[point].area_u
        elseif( grid_region[point].tmask <= 0 ) then
           sea_surface_now[point].u = sea_surface_now[point + {1,0}].t
        elseif( grid_region[point + {1,0}].tmask <= 0) then
           sea_surface_now[point].u = sea_surface_now[point].t
        end


  
    end
  end


end


--------
-------- This writes to a sub field, 2 to N-1, 2 to M
-------- This reads from the full field
--This is the ELEVENTH loop
task update_v_height_launcher( sea_surface_now : region(ispace(int2d), uvt_field), grid_region : region(ispace(int2d), grid_fields) )
                     where writes(sea_surface_now.v), reads(grid_region.area_v, grid_region.area_t, grid_region.tmask, sea_surface_now.t) do


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


  for point in sea_surface_now do
    if( grid_region[point].tmask + grid_region[point+{0,1}].tmask > 0) then
--         if(sshn_t%grid%tmask(ji,jj) * sshn_t%grid%tmask(ji,jj+1) > 0) then
--            rtmp1 = sshn_t%grid%area_t(ji,jj)*sshn_t%data(ji,jj) + &
--                 sshn_t%grid%area_t(ji,jj+1) * sshn_t%data(ji,jj+1)
--            sshn_v%data(ji,jj) = 0.5_go_wp * rtmp1 / sshn_t%grid%area_v(ji,jj)
--         else if(sshn_t%grid%tmask(ji,jj) <= 0) then
--            sshn_v%data(ji,jj) = sshn_t%data(ji,jj+1)
--         else if(sshn_t%grid%tmask(ji,jj+1) <= 0) then
--            sshn_v%data(ji,jj) = sshn_t%data(ji,jj)
--         end if
         if(grid_region[point].tmask * grid_region[point + {0,1}].tmask > 0) then
           var rtmp1 : float = grid_region[point].area_t * sea_surface_now[point].t
                          + grid_region[point+{0,1}].area_t * sea_surface_now[point+{0,1}.t
           sea_surface_now[point].v = 0.5 * rtmp1 / grid_region[point].area_v
        elseif( grid_region[point].tmask <= 0 ) then
           sea_surface_now[point].v = sea_surface_now[point + {0,1}].t
        elseif( grid_region[point + {0,1}].tmask <= 0) then
           sea_surface_now[point].v = sea_surface_now[point].t
        end 

    end

  end



end 
