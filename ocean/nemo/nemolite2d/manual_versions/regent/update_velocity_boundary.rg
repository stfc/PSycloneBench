import "regent"

require("initialise_grid_points")
require("model_init")
local c = regentlib.c

--This is the FIFTH loop
--Loops over 1 to N+1, 1 to M
-- Reads all
task update_uvel_boundary( velocity_after : region(ispace(int2d), uv_field),
                           grid_region : region(ispace(int2d), grid_fields) )
     where writes(velocity_after.u), reads(grid_region.tmask) do
--    do jj = 1, N+1, 1
--       do ji = 1, M, 1
--!          call bc_solid_u_code(ji, jj, &
--!                               ua%data, va%grid%tmask)
--
--!> \todo It's more compiler-friendly to separately compare these two
--!! integer masks with zero but that's a kernel-level optimisation.
--          if(sshn_t%grid%tmask(ji,jj) * sshn_t%grid%tmask(ji+1,jj) == 0)then
--             ua%data(ji,jj) = 0._go_wp
--          end if
--
--       end do
--    end do


  for point in velocity_after do
    if( grid_region[point].tmask * grid_region[point + {1,0}].tmask == int1d(0)) then
      velocity_after[point].u = 0.0
    end
  end

end

--This is the SIXTH loop
--Loops over 1 to N, 1 to M+1
task update_vvel_boundary( velocity_after : region(ispace(int2d), uv_field),
                           grid_region : region(ispace(int2d), grid_fields) )
     where writes(velocity_after.v), reads(grid_region.tmask) do

--    do jj = 1, N, 1
--       do ji = 1, M+1, 1
--!          call bc_solid_u_code(ji, jj, &
--!                               ua%data, va%grid%tmask)
--
--!> \todo It's more compiler-friendly to separately compare these two
--!! integer masks with zero but that's a kernel-level optimisation.
--          if(sshn_t%grid%tmask(ji,jj) * sshn_t%grid%tmask(ji+1,jj) == 0)then
--             va%data(ji,jj) = 0._go_wp
--          end if
--
--       end do
--    end do

  for point in velocity_after do
    if( grid_region[point].tmask * grid_region[point + {0,1}].tmask == int1d(0)) then
      velocity_after[point].v = 0.0
    end
  end


end
