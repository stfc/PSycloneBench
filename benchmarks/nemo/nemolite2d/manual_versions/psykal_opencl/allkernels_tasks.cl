// File originally generated by PSyclone but modified with some manual optimisations:
// - Each kernel iterates over the whole buffer. They are expected to be called with
//   EnqueueTask (or a EnqueueNDRange of just 1 element) instead of a ranged NDRange.
// - LEN is given at compile-time.
// - Added vec_type_hint(double) and xcl_zero_global_work_offset kernels hints.
// - Added const attributes.
// - Removed duplicated LEN expressions.


// LEN value predefined for a problem size of 250x250
#define LEN 256

__attribute__((vec_type_hint(double)))
__attribute__ ((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void bc_flather_u_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict ua,
  const __global double * restrict hu,
  const __global double * restrict sshn_u,
  const __global int * restrict tmask,
  double g
  ){

  // Create buffers for burst copies
  double ua_buffer[LEN];
  double hu_buffer[LEN];
  double sshn_u_buffer[LEN];
  int tmask_buffer[LEN];

  for (int jj = ystart; jj <= ystop; jj++){

      // Burst data reads
      for (int ji = xstart; ji <= xstop; ji++){
        ua_buffer[ji] = ua[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hu_buffer[ji] = hu[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_u_buffer[ji] = sshn_u[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_buffer[ji] = tmask[jj * LEN + ji];
      }

      // Computation loop
      for (int ji = xstart; ji <= xstop; ji++){
          if (!((tmask_buffer[ji] + tmask_buffer[ji + 1]) <= (-1))) {
              if ((tmask[jj * LEN + ji] < 0)) {
                ua_buffer[ji] = (ua_buffer[ji + 1] + (sqrt((g / hu_buffer[ji])) * (sshn_u_buffer[ji] - sshn_u_buffer[ji + 1])));
              } else {
                if ((tmask_buffer[ji + 1] < 0)) {
                  ua_buffer[ji] = (ua_buffer[ji - 1] + (sqrt((g / hu_buffer[ji])) * (sshn_u_buffer[ji] - sshn_u_buffer[ji - 1])));
                }
              }
          }
      }

      // Burst data write
      for (int ji = xstart; ji <= xstop; ji++){
        ua[jj * LEN + ji] = ua_buffer[ji];
      }
  }
}

__attribute__((vec_type_hint(double)))
__attribute__ ((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void bc_flather_v_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict va,
  const __global double * restrict hv,
  const __global double * restrict sshn_v,
  const __global int * restrict tmask,
  double g
  ){
  // Create buffers for burst copies
  double va_buffer[LEN];
  double hv_buffer[LEN];
  double sshn_v_buffer[LEN];
  int tmask_buffer[LEN];

  int tmask_nextrow[LEN];
  double va_previousrow[LEN];
  double va_nextrow[LEN];
  double sshn_v_previousrow[LEN];
  double sshn_v_nextrow[LEN];

  for (int jj = ystart; jj <= ystop; jj++){

      // Burst data reads
      for (int ji = xstart; ji <= xstop; ji++){
        va_buffer[ji] = va[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hv_buffer[ji] = hv[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_v_buffer[ji] = sshn_v[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_buffer[ji] = tmask[jj * LEN + ji];
      }

      // Others
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_nextrow[ji] = tmask[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        va_previousrow[ji] = va[(jj-1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        va_nextrow[ji] = va[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_v_previousrow[ji] = sshn_v[(jj-1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_v_nextrow[ji] = sshn_v[(jj+1) * LEN + ji];
      }

      // Computation loop
      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask_buffer[ji] + tmask_nextrow[ji]) <= (-1))) {
            continue;
          }
          if ((tmask_buffer[ji] < 0)) {
            va_buffer[ji] = (va_nextrow[ji] + (sqrt((g / hv_buffer[ji])) * (sshn_v_buffer[ji] - sshn_v_nextrow[ji])));
          } else {
            if ((tmask_nextrow[ji] < 0)) {
              va_buffer[ji] = (va_previousrow[ji] + (sqrt((g / hv_buffer[ji])) * (sshn_v_buffer[ji] - sshn_v_previousrow[ji])));
            }
          }
      }

      // Burst data writes
      for (int ji = xstart; ji <= xstop; ji++){
        va[jj * LEN + ji] = va_buffer[ji];
      }
  }
}

__attribute__((vec_type_hint(double)))
__attribute__((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void bc_solid_u_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict ua,
  const __global int * restrict tmask
  ){
  for (int jj = ystart; jj <= ystop; jj++){
      double this_tmask;
      double next_tmask = tmask[jj * LEN + xstart];
      for (int ji = xstart; ji <= xstop; ji++){
          this_tmask = next_tmask;
          next_tmask = tmask[jj * LEN + (ji + 1)];
          if (((this_tmask * next_tmask) == 0)) {
            ua[jj * LEN + ji] = 0.;
          }
      }
  }
}

__attribute__((vec_type_hint(double)))
__attribute__ ((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void bc_solid_v_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict va,
  const __global int * restrict tmask
  ){
  double this_row[LEN];
  double next_row[LEN];
  for (int ji = xstart; ji <= xstop; ji++){
      next_row[ji] = tmask[ystart * LEN + ji];
  }

  for (int jj = ystart; jj <= ystop; jj++){
      for (int ji = xstart; ji <= xstop; ji++){
          this_row[ji] = next_row[ji];
          next_row[ji] = tmask[(jj + 1) * LEN + ji];
      }

      for (int ji = xstart; ji <= xstop; ji++){
          if ((this_row[ji] * next_row[jj]) == 0) {
            va[jj * LEN + ji] = 0.;
          }
      }
  }
}

__attribute__((vec_type_hint(double)))
__attribute__ ((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void bc_ssh_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  int istep,
  __global double * restrict ssha,
  const __global int * restrict tmask,
  double rdt
  ){

  double amp_tide = 0.2;
  double omega_tide = ((2.0 * 3.14159) / (12.42 * 3600.));
  double rtime = ((double)istep * rdt);
  double value = (amp_tide * sin((omega_tide * rtime)));

  // Temporal buffers
  int tmask_p1[LEN];
  int tmask_this[LEN];
  int tmask_m1[LEN];
  
  // Next: can also buffer ssha
  for (int ji = xstart; ji <= xstop; ji++){
    tmask_m1[ji] = tmask[(ystart - 1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    tmask_this[ji] = tmask[ystart * LEN + ji];
  }


  for (int jj = ystart; jj <= ystop; jj++){

      for (int ji = xstart; ji <= xstop; ji++){
        tmask_p1[ji] = tmask[(jj+1) * LEN + ji];
      }

      for (int ji = xstart; ji <= xstop; ji++){
          if (!(tmask[jj * LEN + ji] <= 0)) {
              if ((tmask_m1[ji] < 0)) {
                ssha[jj * LEN + ji] = value;
              } else {
                if ((tmask_p1[ji] < 0)) {
                  ssha[jj * LEN + ji] = value;
                } else {
                  if ((tmask_this[ji+1] < 0)) {
                    ssha[jj * LEN + ji] = value;
                  } else {
                    if ((tmask_this[ji-1] < 0)) {
                      ssha[jj * LEN + ji] = value;
                    }
                  }
                }
              }
          }
      }

      for (int ji = xstart; ji <= xstop; ji++){
        tmask_m1[ji] = tmask_this[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_this[ji] = tmask_p1[ji];
      }

  }
}

__attribute__((vec_type_hint(double)))
__attribute__ ((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void continuity_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict ssha,
  const __global double * restrict sshn,
  const __global double * restrict sshn_u,
  const __global double * restrict sshn_v,
  const __global double * restrict hu,
  const __global double * restrict hv,
  const __global double * restrict un,
  const __global double * restrict vn,
  const __global double * restrict e12t,
  double rdt
  ){
  // Create buffers for burst copies
  double ssha_buffer[LEN];
  double sshn_buffer[LEN];
  double e12t_buffer[LEN];

  double sshn_u_buffer[LEN];
  double hu_buffer[LEN];
  double un_buffer[LEN];
  double sshn_v_buffer[LEN];
  double hv_buffer[LEN];
  double vn_buffer[LEN];

  // Create row cache to to save a row for the next jj iteration
  double ssh_v_previousrow[LEN];
  double hv_previousrow[LEN];
  double vn_previousrow[LEN];

  // Burst load previous rows
  for (int ji = xstart; ji <= xstop; ji++){
    ssh_v_previousrow[ji] = sshn_v[(ystart - 1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    hv_previousrow[ji] = hv[(ystart - 1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    vn_previousrow[ji] = vn[(ystart - 1) * LEN + ji];
  }

  for (int jj = ystart; jj <= ystop; jj++){
      // Keep a copy of this 3 values for the next ji iteration
      double sshn_u_jim1 = sshn_u[jj * LEN + (xstart - 1)];
      double hu_jim1 = hu[jj * LEN + (xstart - 1)];
      double un_jim1 = un[jj * LEN + (xstart - 1)];

      for (int ji = xstart; ji <= xstop; ji++){
        sshn_buffer[ji] = sshn[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e12t_buffer[ji] = e12t[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_u_buffer[ji] = sshn_u[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hu_buffer[ji] = hu[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        un_buffer[ji] = un[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_v_buffer[ji] = sshn_v[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hv_buffer[ji] = hv[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        vn_buffer[ji] = vn[jj * LEN + ji];
      }

      for (int ji = xstart; ji <= xstop; ji++){
          double this_sshn_u = sshn_u_buffer[ji];
          double this_hu = hu_buffer[ji];
          double this_un = un_buffer[ji];

          double this_sshn_v = sshn_v_buffer[ji];
          double this_hv = hv_buffer[ji];
          double this_vn = vn_buffer[ji];

          // Kernel computation
          double rtmp1 = ((this_sshn_u + this_hu) * this_un);
          double rtmp2 = ((sshn_u_jim1 + hu_jim1) * un_jim1);
          double rtmp3 = ((this_sshn_v + this_hv) * this_vn);
          double rtmp4 = ((ssh_v_previousrow[ji] + hv_previousrow[ji]) * vn_previousrow[ji]);
          ssha_buffer[ji] = (sshn_buffer[ji] + (((((rtmp2 - rtmp1) + rtmp4) - rtmp3) * rdt) / e12t_buffer[ji]));

          // Save the current values in registers for next iteration
          sshn_u_jim1 = this_sshn_u;
          hu_jim1 = this_hu;
          un_jim1 = this_un;
          ssh_v_previousrow[ji] = this_sshn_v; 
          hv_previousrow[ji] = this_hv;
          vn_previousrow[ji] = this_vn;
      }

      for (int ji = xstart; ji <= xstop; ji++){
        ssha[jj * LEN + ji] = ssha_buffer[ji];
      }
  }
}

__attribute__((vec_type_hint(double)))
__attribute__ ((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void field_copy_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict output,
  const __global double * restrict input
  ){
  for (int jj = ystart; jj <= ystop; jj++){
      for (int ji = xstart; ji <= xstop; ji++){
        output[jj * LEN + ji] = input[jj * LEN + ji];
      }
  }
}

__attribute__((vec_type_hint(double)))
__attribute__ ((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void momentum_u_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict ua,
  const __global double * restrict un,
  const __global double * restrict vn,
  const __global double * restrict hu,
  const __global double * restrict hv,
  const __global double * restrict ht,
  const __global double * restrict ssha_u,
  const __global double * restrict sshn,
  const __global double * restrict sshn_u,
  const __global double * restrict sshn_v,
  const __global int * restrict tmask,
  const __global double * restrict e1u,
  const __global double * restrict e1v,
  const __global double * restrict e1t,
  const __global double * restrict e2u,
  const __global double * restrict e2t,
  const __global double * restrict e12u,
  const __global double * restrict gphiu,
  double omega,
  double d2r,
  double g,
  double rdt,
  double cbfr,
  double visc
  ){
  double u_e;
  double u_w;
  double v_n;
  double v_s;
  double v_nc;
  double v_sc;
  double depe;
  double depw;
  double deps;
  double depn;
  double hpg;
  double adv;
  double cor;
  double vis;
  double dudx_e;
  double dudx_w;
  double dudy_s;
  double dudy_n;
  double uu_e;
  double uu_n;
  double uu_s;
  double uu_w;

  // Create buffers for burst copies
  double ua_buffer[LEN];
  double un_buffer[LEN];
  //#pragma HLS array_partition variable=un_buffer cyclic factor=3
  double vn_buffer[LEN];
  double hu_buffer[LEN];
  double hv_buffer[LEN];
  double ht_buffer[LEN];
  double ssha_u_buffer[LEN];
  double sshn_buffer[LEN];
  double sshn_u_buffer[LEN];
  double sshn_v_buffer[LEN];
  int tmask_buffer[LEN];
  double e1u_buffer[LEN];
  double e1v_buffer[LEN];
  double e1t_buffer[LEN];
  double e2u_buffer[LEN];
  double e2t_buffer[LEN];
  double e12u_buffer[LEN];
  double gphiu_buffer[LEN];

  double un_previousrow[LEN];
  double un_nextrow[LEN];
  double vn_previousrow[LEN];
  double hu_previousrow[LEN];
  double hu_nextrow[LEN];
  double hv_previousrow[LEN];
  double sshn_u_previousrow[LEN];
  double sshn_u_nextrow[LEN];
  double sshn_v_previousrow[LEN];
  int tmask_previousrow[LEN];
  int tmask_nextrow[LEN];
  int e1v_previousrow[LEN];
  int e2u_previousrow[LEN];
  int e2u_nextrow[LEN];

  // Load initial data for nextrow that are needed to populate buffers
  for (int ji = xstart; ji <= xstop; ji++){
    un_nextrow[ji] = un[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    hu_nextrow[ji] = hu[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    tmask_nextrow[ji] = tmask[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    e2u_nextrow[ji] = e2u[ystart * LEN + ji];
  }
  // Load initial data for buffers that are needed to populate previousrows
  for (int ji = xstart; ji <= xstop; ji++){
    un_buffer[ji] = un[(ystart-1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    vn_buffer[ji] = vn[(ystart-1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    hu_buffer[ji] = hu[(ystart-1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    hv_buffer[ji] = hv[(ystart-1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    sshn_u_buffer[ji] = sshn_u[(ystart-1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    sshn_v_buffer[ji] = sshn_v[(ystart-1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    tmask_buffer[ji] = tmask[(ystart-1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    e1v_buffer[ji] = e1v[(ystart-1) * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    e2u_buffer[ji] = e2u[(ystart-1) * LEN + ji];
  }

  for (int jj = ystart; jj <= ystop; jj++){

      // Populate previous rows (from already loaded data in buffers if possible)
      for (int ji = xstart; ji <= xstop; ji++){
        un_previousrow[ji] = un_buffer[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        vn_previousrow[ji] = vn_buffer[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hu_previousrow[ji] = hu_buffer[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hv_previousrow[ji] = hv_buffer[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_u_previousrow[ji] = sshn_u_buffer[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_u_nextrow[ji] = sshn_u_buffer[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_v_previousrow[ji] = sshn_v_buffer[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_previousrow[ji] = tmask_buffer[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e1v_previousrow[ji] = e1v_buffer[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e2u_previousrow[ji] = e2u_buffer[ji];
      }

      // Burst data reads (from already loaded data in nextrows if possible)
      for (int ji = xstart; ji <= xstop; ji++){
        ua_buffer[ji] = ua[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        un_buffer[ji] = un_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        vn_buffer[ji] = vn[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hu_buffer[ji] = hu_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hv_buffer[ji] = hv[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        ht_buffer[ji] = ht[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        ssha_u_buffer[ji] = ssha_u[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_buffer[ji] = sshn[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_u_buffer[ji] = sshn_u[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_v_buffer[ji] = sshn_v[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_buffer[ji] = tmask_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e1u_buffer[ji] = e1u[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e1v_buffer[ji] = e1v[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e1t_buffer[ji] = e1t[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e2u_buffer[ji] = e2u_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e2t_buffer[ji] = e2t[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e12u_buffer[ji] = e12u[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        gphiu_buffer[ji] = gphiu[jj * LEN + ji];
      }

      // Load next row with data burst reads
      for (int ji = xstart; ji <= xstop; ji++){
        un_nextrow[ji] = un[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hu_nextrow[ji] = hu[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_nextrow[ji] = tmask[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e2u_nextrow[ji] = e2u[(jj+1) * LEN + ji];
      }

      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask_buffer[ji] + tmask_buffer[(ji + 1)]) <= 0)) {
            continue;
          }
          if (((tmask_buffer[ji] <= 0) || (tmask_buffer[(ji + 1)] <= 0))) {
            continue;
          }
          u_e = ((0.5 * (un_buffer[ji] + un_buffer[(ji + 1)])) * e2t_buffer[(ji + 1)]);
          depe = (ht_buffer[(ji + 1)] + sshn_buffer[(ji + 1)]);
          u_w = ((0.5 * (un_buffer[ji] + un_buffer[(ji - 1)])) * e2t_buffer[ji]);
          depw = (ht_buffer[ji] + sshn_buffer[ji]);
          v_sc = (0.5 * (vn_previousrow[ji] + vn_previousrow[(ji + 1)]));
          v_s = ((0.5 * v_sc) * (e1v_previousrow[ji] + e1v_previousrow[(ji + 1)]));
          deps = (0.5 * (((hv_previousrow[ji] + sshn_v_previousrow[ji]) + hv_previousrow[(ji + 1)]) + sshn_v_previousrow[(ji + 1)]));
          v_nc = (0.5 * (vn_buffer[ji] + vn_buffer[(ji + 1)]));
          v_n = ((0.5 * v_nc) * (e1v_buffer[ji] + e1v_buffer[(ji + 1)]));
          depn = (0.5 * (((hv_buffer[ji] + sshn_v_buffer[ji]) + hv_buffer[(ji + 1)]) + sshn_v_buffer[(ji + 1)]));
          uu_w = (((0.5 - copysign(0.5, u_w)) * un_buffer[ji]) + ((0.5 + copysign(0.5, u_w)) * un_buffer[(ji - 1)]));
          uu_e = (((0.5 + copysign(0.5, u_e)) * un_buffer[ji]) + ((0.5 - copysign(0.5, u_e)) * un_buffer[(ji + 1)]));
          if (((tmask_previousrow[ji] <= 0) || (tmask_previousrow[(ji + 1)] <= 0))) {
            uu_s = ((0.5 - copysign(0.5, v_s)) * un_buffer[ji]);
          } else {
            uu_s = (((0.5 - copysign(0.5, v_s)) * un_buffer[ji]) + ((0.5 + copysign(0.5, v_s)) * un_previousrow[ji]));
          }
          if (((tmask_nextrow[ji] <= 0) || (tmask_nextrow[(ji + 1)] <= 0))) {
            uu_n = ((0.5 + copysign(0.5, v_n)) * un_buffer[ji]);
          } else {
            uu_n = (((0.5 + copysign(0.5, v_n)) * un_buffer[ji]) + ((0.5 - copysign(0.5, v_n)) * un_nextrow[ji]));
          }
          adv = (((((uu_w * u_w) * depw) - ((uu_e * u_e) * depe)) + ((uu_s * v_s) * deps)) - ((uu_n * v_n) * depn));
          dudx_e = (((un_buffer[(ji + 1)] - un_buffer[ji]) / e1t_buffer[(ji + 1)]) * (ht_buffer[(ji + 1)] + sshn_buffer[(ji + 1)]));
          dudx_w = (((un_buffer[ji] - un_buffer[(ji - 1)]) / e1t_buffer[ji]) * (ht_buffer[ji] + sshn_buffer[ji]));
          if (((tmask_previousrow[ji] <= 0) || (tmask_previousrow[(ji + 1)] <= 0))) {
            dudy_s = 0.0;
          } else {
            dudy_s = (((un_buffer[ji] - un_previousrow[ji]) / (e2u_buffer[ji] + e2u_previousrow[ji])) * (((hu_buffer[ji] + sshn_u_buffer[ji]) + hu_previousrow[ji]) + sshn_u_previousrow[ji]));
          }
          if (((tmask_nextrow[ji] <= 0) || (tmask_nextrow[(ji + 1)] <= 0))) {
            dudy_n = 0.0;
          } else {
            dudy_n = (((un_nextrow[ji] - un_buffer[ji]) / (e2u_buffer[ji] + e2u_nextrow[ji])) * (((hu_buffer[ji] + sshn_u_buffer[ji]) + hu_nextrow[ji]) + sshn_u_nextrow[ji]));
          }
          vis = (((dudx_e - dudx_w) * e2u_buffer[ji]) + (((dudy_n - dudy_s) * e1u_buffer[ji]) * 0.5));
          vis = (visc * vis);
          cor = (((0.5 * (((2. * omega) * sin((gphiu_buffer[ji] * d2r))) * (v_sc + v_nc))) * e12u_buffer[ji]) * (hu_buffer[ji] + sshn_u_buffer[ji]));
          hpg = (-(((g * (hu_buffer[ji] + sshn_u_buffer[ji])) * e2u_buffer[ji]) * (sshn_buffer[(ji + 1)] - sshn_buffer[ji])));
          ua_buffer[ji] = ((((un_buffer[ji] * (hu_buffer[ji] + sshn_u_buffer[ji])) + ((rdt * (((adv + vis) + cor) + hpg)) / e12u_buffer[ji])) / (hu_buffer[ji] + ssha_u_buffer[ji])) / (1.0 + (cbfr * rdt)));
      }
      for (int ji = xstart; ji <= xstop; ji++){
        ua[jj * LEN + ji] = ua_buffer[ji];
      }
  }
}

__attribute__((vec_type_hint(double)))
__attribute__ ((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void momentum_v_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict va,
  const __global double * restrict un,
  const __global double * restrict vn,
  const __global double * restrict hu,
  const __global double * restrict hv,
  const __global double * restrict ht,
  const __global double * restrict ssha_v,
  const __global double * restrict sshn,
  const __global double * restrict sshn_u,
  const __global double * restrict sshn_v,
  const __global int * restrict tmask,
  const __global double * restrict e1v,
  const __global double * restrict e1t,
  const __global double * restrict e2u,
  const __global double * restrict e2v,
  const __global double * restrict e2t,
  const __global double * restrict e12v,
  const __global double * restrict gphiv,
  double omega,
  double d2r,
  double g,
  double rdt,
  double cbfr,
  double visc
  ){
  double u_e;
  double u_w;
  double v_n;
  double v_s;
  double u_ec;
  double u_wc;
  double vv_e;
  double vv_n;
  double vv_s;
  double vv_w;
  double depe;
  double depw;
  double deps;
  double depn;
  double hpg;
  double adv;
  double cor;
  double vis;
  double dvdx_e;
  double dvdx_w;
  double dvdy_n;
  double dvdy_s;

  // Create buffers for burst copies
  double va_buffer[LEN];
  double un_buffer[LEN];
  double vn_buffer[LEN];
  double hu_buffer[LEN];
  double hv_buffer[LEN];
  double ht_buffer[LEN];
  double ssha_v_buffer[LEN];
  double sshn_buffer[LEN];
  double sshn_u_buffer[LEN];
  double sshn_v_buffer[LEN];
  int tmask_buffer[LEN];
  double e1v_buffer[LEN];
  double e1t_buffer[LEN];
  double e2u_buffer[LEN];
  double e2v_buffer[LEN];
  double e2t_buffer[LEN];
  double e12v_buffer[LEN];
  double gphiv_buffer[LEN];

  double un_nextrow[LEN];
  double vn_previousrow[LEN];
  double vn_nextrow[LEN];
  double hu_nextrow[LEN];
  double ht_nextrow[LEN];
  double sshn_nextrow[LEN];
  double sshn_u_nextrow[LEN];
  int tmask_nextrow[LEN];
  double e1t_nextrow[LEN];
  double e2u_nextrow[LEN];
  double e2t_nextrow[LEN];

  // Load initial data for nextrow that are needed to populate buffers
  for (int ji = xstart; ji <= xstop; ji++){
    un_nextrow[ji] = un[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    vn_nextrow[ji] = vn[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    hu_nextrow[ji] = hu[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    ht_nextrow[ji] = ht[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    sshn_nextrow[ji] = sshn[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    sshn_u_nextrow[ji] = sshn_u[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    tmask_nextrow[ji] = tmask[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    e1t_nextrow[ji] = e1t[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    e2u_nextrow[ji] = e2u[ystart * LEN + ji];
  }
  for (int ji = xstart; ji <= xstop; ji++){
    e2t_nextrow[ji] = e2t[ystart * LEN + ji];
  }

  // Load initial data for buffers that are needed to populate previousrows
  for (int ji = xstart; ji <= xstop; ji++){
    vn_buffer[ji] = vn[(ystart-1) * LEN + ji];
  }
  
  for (int jj = ystart; jj <= ystop; jj++){

      for (int ji = xstart; ji <= xstop; ji++){
        vn_previousrow[ji] = vn_buffer[ji];
      }

      // Burst buffer reads
      for (int ji = xstart; ji <= xstop; ji++){
        va_buffer[ji] = va[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        un_buffer[ji] = un_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        vn_buffer[ji] = vn_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hu_buffer[ji] = hu_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hv_buffer[ji] = hv[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        ht_buffer[ji] = ht_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        ssha_v_buffer[ji] = ssha_v[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_buffer[ji] = sshn_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_u_buffer[ji] = sshn_u_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_v_buffer[ji] = sshn_v[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_buffer[ji] = tmask_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e1v_buffer[ji] = e1v[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e1t_buffer[ji] = e1t_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e2u_buffer[ji] = e2u_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e2v_buffer[ji] = e2v[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e2t_buffer[ji] = e2t_nextrow[ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e12v_buffer[ji] = e12v[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        gphiv_buffer[ji] = gphiv[jj * LEN + ji];
      }


      for (int ji = xstart; ji <= xstop; ji++){
        un_nextrow[ji] = un[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        vn_nextrow[ji] = vn[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        hu_nextrow[ji] = hu[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        ht_nextrow[ji] = ht[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_nextrow[ji] = sshn[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_u_nextrow[ji] = sshn_u[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_nextrow[ji] = tmask[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e1t_nextrow[ji] = e1t[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e2u_nextrow[ji] = e2u[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e2t_nextrow[ji] = e2t[(jj+1) * LEN + ji];
      }

      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask_buffer[ji] + tmask_buffer[(ji + 1)]) <= 0)) {
            continue;
          }
          if (((tmask_buffer[ji] <= 0) || (tmask_nextrow[ji] <= 0))) {
            continue;
          }
          v_n = ((0.5 * (vn_buffer[ji] + vn_nextrow[ji])) * e1t_nextrow[ji]);
          depn = (ht_nextrow[ji] + sshn_nextrow[ji]);
          v_s = ((0.5 * (vn_buffer[ji] + vn_previousrow[ji])) * e1t_buffer[ji]);
          deps = (ht_nextrow[ji] + sshn_buffer[ji]);
          u_wc = (0.5 * (un_buffer[(ji - 1)] + un_nextrow[(ji - 1)]));
          u_w = ((0.5 * u_wc) * (e2u_buffer[(ji - 1)] + e2u_nextrow[(ji - 1)]));
          depw = (0.50 * (((hu_buffer[(ji - 1)] + sshn_u_buffer[(ji - 1)]) + hu_nextrow[(ji - 1)]) + sshn_u_nextrow[(ji - 1)]));
          u_ec = (0.5 * (un_buffer[ji] + un_nextrow[ji]));
          u_e = ((0.5 * u_ec) * (e2u_buffer[ji] + e2u_nextrow[ji]));
          depe = (0.50 * (((hu_buffer[ji] + sshn_u_buffer[ji]) + hu_nextrow[ji]) + sshn_u_nextrow[ji]));
          vv_s = (((0.5 - copysign(0.5, v_s)) * vn_buffer[ji]) + ((0.5 + copysign(0.5, v_s)) * vn_previousrow[ji]));
          vv_n = (((0.5 + copysign(0.5, v_n)) * vn_buffer[ji]) + ((0.5 - copysign(0.5, v_n)) * vn_nextrow[ji]));
          if (((tmask_buffer[(ji - 1)] <= 0) || (tmask_nextrow[(ji - 1)] <= 0))) {
            vv_w = ((0.5 - copysign(0.5, u_w)) * vn_buffer[ji]);
          } else {
            vv_w = (((0.5 - copysign(0.5, u_w)) * vn_buffer[ji]) + ((0.5 + copysign(0.5, u_w)) * vn_buffer[(ji - 1)]));
          }
          if (((tmask_buffer[(ji + 1)] <= 0) || (tmask_nextrow[(ji + 1)] <= 0))) {
            vv_e = ((0.5 + copysign(0.5, u_e)) * vn_buffer[ji]);
          } else {
            vv_e = (((0.5 + copysign(0.5, u_e)) * vn_buffer[ji]) + ((0.5 - copysign(0.5, u_e)) * vn_buffer[(ji + 1)]));
          }
          adv = (((((vv_w * u_w) * depw) - ((vv_e * u_e) * depe)) + ((vv_s * v_s) * deps)) - ((vv_n * v_n) * depn));
          dvdy_n = (((vn_nextrow[ji] - vn_buffer[ji]) / e2t_nextrow[ji]) * (ht_nextrow[ji] + sshn_nextrow[ji]));
          dvdy_s = (((vn_buffer[ji] - vn_previousrow[ji]) / e2t_buffer[ji]) * (ht_buffer[ji] + sshn_buffer[ji]));
          if (((tmask_buffer[(ji - 1)] <= 0) || (tmask_nextrow[(ji - 1)] <= 0))) {
            dvdx_w = 0.0;
          } else {
            dvdx_w = (((vn_buffer[ji] - vn_buffer[(ji - 1)]) / (e1v_buffer[ji] + e1v_buffer[(ji - 1)])) * (((hv_buffer[ji] + sshn_v_buffer[ji]) + hv_buffer[(ji - 1)]) + sshn_v_buffer[(ji - 1)]));
          }
          if (((tmask_buffer[(ji + 1)] <= 0) || (tmask_nextrow[(ji + 1)] <= 0))) {
            dvdx_e = 0.0;
          } else {
            dvdx_e = (((vn_buffer[(ji + 1)] - vn_buffer[ji]) / (e1v_buffer[ji] + e1v_buffer[(ji + 1)])) * (((hv_buffer[ji] + sshn_v_buffer[ji]) + hv_buffer[(ji + 1)]) + sshn_v_buffer[(ji + 1)]));
          }
          vis = (((dvdy_n - dvdy_s) * e1v_buffer[ji]) + (((dvdx_e - dvdx_w) * e2v_buffer[ji]) * 0.5));
          vis = (visc * vis);
          cor = (-(((0.5 * (((2. * omega) * sin((gphiv_buffer[ji] * d2r))) * (u_ec + u_wc))) * e12v_buffer[ji]) * (hv_buffer[ji] + sshn_v_buffer[ji])));
          hpg = (-(((g * (hv_buffer[ji] + sshn_v_buffer[ji])) * e1v_buffer[ji]) * (sshn_nextrow[ji] - sshn_buffer[ji])));
          va_buffer[ji] = ((((vn_buffer[ji] * (hv_buffer[ji] + sshn_v_buffer[ji])) + ((rdt * (((adv + vis) + cor) + hpg)) / e12v_buffer[ji])) / (hv_buffer[ji] + ssha_v_buffer[ji])) / (1.0 + (cbfr * rdt)));
      }

      for (int ji = xstart; ji <= xstop; ji++){
        va[jj * LEN + ji] = va_buffer[ji];
      }

  }
}

__attribute__((vec_type_hint(double)))
__attribute__ ((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void next_sshu_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict sshn_u,
  const __global double * restrict sshn,
  const __global int * restrict tmask,
  const __global double * restrict e12t,
  const __global double * restrict e12u
  ){
  // Create buffers for burst copies
  double sshn_u_buffer[LEN];
  double sshn_buffer[LEN];
  double e12t_buffer[LEN];
  double e12u_buffer[LEN];
  int tmask_buffer[LEN];

  for (int jj = ystart; jj <= ystop; jj++){

      // Burst data reads
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_buffer[ji] = sshn[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_buffer[ji] = tmask[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e12t_buffer[ji] = e12t[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e12u_buffer[ji] = e12u[jj * LEN + ji];
      }

      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask_buffer[ji] + tmask_buffer[(ji + 1)]) <= 0)) {
            continue;
          }
          if (((tmask_buffer[ji] * tmask_buffer[(ji + 1)]) > 0)) {
            double rtmp1 = ((e12t_buffer[ji] * sshn_buffer[ji]) + (e12t_buffer[(ji + 1)] * sshn_buffer[(ji + 1)]));
            sshn_u_buffer[ji] = ((0.5 * rtmp1) / e12u_buffer[ji]);
          } else {
            if ((tmask_buffer[ji] <= 0)) {
              sshn_u_buffer[ji] = sshn_buffer[(ji + 1)];
            } else {
              if ((tmask[jj * LEN + (ji + 1)] <= 0)) {
                sshn_u_buffer[ji] = sshn_buffer[ji];
              }
            }
          }
      }
      // Burst data write
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_u[jj * LEN + ji] = sshn_u_buffer[ji];
      }
  }
}

__attribute__((vec_type_hint(double)))
__attribute__ ((reqd_work_group_size(1, 1, 1)))
__attribute__((xcl_zero_global_work_offset))
__kernel void next_sshv_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict sshn_v,
  const __global double * restrict sshn,
  const __global int * restrict tmask,
  const __global double * restrict e12t,
  const __global double * restrict e12v
  ){
  // Create buffers for burst copies
  double sshn_v_buffer[LEN];
  double sshn_buffer[LEN];
  double e12t_buffer[LEN];
  double e12v_buffer[LEN];
  int tmask_buffer[LEN];

  double sshn_nextrow[LEN];
  int tmask_nextrow[LEN];
  double e12t_nextrow[LEN];

  for (int jj = ystart; jj <= ystop; jj++){

      // Burst data reads
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_buffer[ji] = sshn[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_buffer[ji] = tmask[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e12t_buffer[ji] = e12t[jj * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e12v_buffer[ji] = e12v[jj * LEN + ji];
      }

      for (int ji = xstart; ji <= xstop; ji++){
        sshn_nextrow[ji] = sshn[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        tmask_nextrow[ji] = tmask[(jj+1) * LEN + ji];
      }
      for (int ji = xstart; ji <= xstop; ji++){
        e12t_nextrow[ji] = e12t[(jj+1) * LEN + ji];
      }
        

      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask_buffer[ji] + tmask_nextrow[ji]) <= 0)) {
            continue;
          }
          if (((tmask_buffer[ji] * tmask_nextrow[ji]) > 0)) {
            double rtmp1 = ((e12t_buffer[ji] * sshn_buffer[ji]) + (e12t_nextrow[ji] * sshn_nextrow[ji]));
            sshn_v_buffer[ji] = ((0.5 * rtmp1) / e12v_buffer[ji]);
          } else {
            if ((tmask_buffer[ji] <= 0)) {
              sshn_v_buffer[ji] = sshn_nextrow[ji];
            } else {
              if ((tmask_nextrow[ji] <= 0)) {
                sshn_v_buffer[ji] = sshn_buffer[ji];
              }
            }
          }
      }

      // Burst data write
      for (int ji = xstart; ji <= xstop; ji++){
        sshn_v[jj * LEN + ji] = sshn_v_buffer[ji];
      }
  }
}

