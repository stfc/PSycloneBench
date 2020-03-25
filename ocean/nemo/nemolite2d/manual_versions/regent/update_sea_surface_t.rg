import "regent"


require("model_init")
require("initialise_grid_points")

local c = regentlib.c

sin = c.sin

--This is the FOURTH loop
--Writes to 2 to N 2 to M
-- Reads from all
task update_sea_surface_t(sea_surface_after : region(ispace(int2d),uvt_field),
                          grid_region : region(ispace(int2d), grid_fields),
                          rdt : double, step : int32 )
     where writes(sea_surface_after.t), reads(grid_region.tmask) do


--    DO jj = 2, N
--! SIMD
--       DO ji = 2, M
--!          call bc_ssh_code(ji, jj, &
--!                           istp, ssha%data, sshn_t%grid%tmask)
--
--          amp_tide   = 0.2_go_wp
--          omega_tide = 2.0_go_wp * 3.14159_go_wp / (12.42_go_wp * 3600._go_wp)
--          rtime = real(istp, go_wp) * rdt
--
--          if(sshn_t%grid%tmask(ji,jj) <= 0) cycle
--
--          IF     (sshn_t%grid%tmask(ji,jj-1) < 0) THEN
--             ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
--          ELSE IF(sshn_t%grid%tmask(ji,jj+1) < 0) THEN
--             ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
--          ELSE IF(sshn_t%grid%tmask(ji+1,jj) < 0) THEN
--             ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
--          ELSE IF(sshn_t%grid%tmask(ji-1,jj) < 0) THEN
--             ssha%data(ji,jj) = amp_tide * sin(omega_tide * rtime)
--          END IF
--
--       END DO
--    END DO
  var fstep : double = double(step)
  var amp_tide : double = 0.2
  var omega_tide : double = 2.0 * 3.14159 / (12.42 * 3600.0)
  var rtime : double = fstep * rdt
  var tide_value : double = amp_tide * sin(omega_tide * rtime)
  for point in sea_surface_after do
    if( grid_region[point + {0,-1}].tmask < int1d(0)
      or grid_region[point + {0,1}].tmask < int1d(0)
      or grid_region[point + {1,0}].tmask < int1d(0)
      or grid_region[point + {-1,0}].tmask < int1d(0)) then
      sea_surface_after[point].t = tide_value
    end
  end


end

