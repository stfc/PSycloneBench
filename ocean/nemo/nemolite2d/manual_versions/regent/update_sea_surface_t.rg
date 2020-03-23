import "regent"


requires("model_init")
requires("initialise_grid_points")

local c = regentlib.c

sin = c.sinf

--This is the FOURTH loop
--Writes to 2 to N 2 to M
-- Reads from all
task update_sea_surface_t(sea_surface_after : region(ispace(int2d),uvt_field),
                          grid_region : region(ispace(int2d), grid_fields),
                          rdt : float, step : float )
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

  var amp_tide : float = 0.2
  var omega_tide : float = 0.2 * 3.14159 / (12.42 * 3600.0)
  var rtime = step * rdt
  var tide_value = amp_tide * sin(omega_tide * rtime)
  for point in sea_surface_next do
    if( grid_region[point + {0,-1}].tmask < 0
      or grid_region[point + {0,1}].tmask < 0
      or grid_region[point + {1,0}].tmask < 0
      or grid_region[point + {-1,0}].tmask < 0) then
      sea_surface_after[point] = tide_value
    end
  end


end

