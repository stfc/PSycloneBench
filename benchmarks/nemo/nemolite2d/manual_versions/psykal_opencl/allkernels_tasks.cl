
#define LEN 256

__attribute__((vec_type_hint(double)))
__kernel __attribute__ ((reqd_work_group_size(1, 1, 1)))
__kernel __attribute__((xcl_zero_global_work_offset))
__kernel void bc_flather_u_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict ua,
  __global double * restrict hu,
  __global double * restrict sshn_u,
  __global int * restrict tmask,
  double g
  ){
  
  for (int jj = ystart; jj <= ystop; jj++){
      __attribute__((opencl_unroll_hint(8)))
      for (int ji = xstart; ji <= xstop; ji++){
          int jiu;
          if (!((tmask[jj * LEN + ji] + tmask[jj * LEN + (ji + 1)]) <= (-1))) {
              if ((tmask[jj * LEN + ji] < 0)) {
                jiu = (ji + 1);
                ua[jj * LEN + ji] = (ua[jj * LEN + jiu] + (native_sqrt((g / hu[jj * LEN + ji])) * (sshn_u[jj * LEN + ji] - sshn_u[jj * LEN + jiu])));
              } else {
                if ((tmask[jj * LEN + (ji + 1)] < 0)) {
                  jiu = (ji - 1);
                  ua[jj * LEN + ji] = (ua[jj * LEN + jiu] + (native_sqrt((g / hu[jj * LEN + ji])) * (sshn_u[jj * LEN + ji] - sshn_u[jj * LEN + jiu])));
                }
              }
          }
      }
  }
}

__attribute__((vec_type_hint(double)))
__kernel __attribute__ ((reqd_work_group_size(1, 1, 1)))
__kernel __attribute__((xcl_zero_global_work_offset))
__kernel void bc_flather_v_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict va,
  __global double * restrict hv,
  __global double * restrict sshn_v,
  __global int * restrict tmask,
  double g
  ){
  for (int jj = ystart; jj <= ystop; jj++){
      __attribute__((xcl_pipeline_loop))
      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask[jj * LEN + ji] + tmask[(jj + 1) * LEN + ji]) <= (-1))) {
            continue;
          }
          int jiv;
          if ((tmask[jj * LEN + ji] < 0)) {
            jiv = (jj + 1);
            va[jj * LEN + ji] = (va[jiv * LEN + ji] + (sqrt((g / hv[jj * LEN + ji])) * (sshn_v[jj * LEN + ji] - sshn_v[jiv * LEN + ji])));
          } else {
            if ((tmask[(jj + 1) * LEN + ji] < 0)) {
              jiv = (jj - 1);
              va[jj * LEN + ji] = (va[jiv * LEN + ji] + (sqrt((g / hv[jj * LEN + ji])) * (sshn_v[jj * LEN + ji] - sshn_v[jiv * LEN + ji])));
            }
          }
      }
  }
}

__kernel void bc_solid_u_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict ua,
  __global int * restrict tmask
  ){
  for (int jj = ystart; jj <= ystop; jj++){
      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask[jj * LEN + ji] * tmask[jj * LEN + (ji + 1)]) == 0)) {
            ua[jj * LEN + ji] = 0.;
          }
      }
  }
}

__kernel void bc_solid_v_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict va,
  __global int * restrict tmask
  ){
  for (int jj = ystart; jj <= ystop; jj++){
      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask[jj * LEN + ji] * tmask[(jj + 1) * LEN + ji]) == 0)) {
            va[jj * LEN + ji] = 0.;
          }
      }
  }
}

__attribute__((vec_type_hint(double)))
__kernel __attribute__ ((reqd_work_group_size(1, 1, 1)))
__kernel __attribute__((xcl_zero_global_work_offset))
__kernel void bc_ssh_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  int istep,
  __global double * restrict ssha,
  __global int * restrict tmask,
  double rdt
  ){
  for (int jj = ystart; jj <= ystop; jj++){
      __attribute__((opencl_unroll_hint(8)))
      for (int ji = xstart; ji <= xstop; ji++){
          double amp_tide = 0.2;
          double omega_tide = ((2.0 * 3.14159) / (12.42 * 3600.));
          double rtime = ((double)istep * rdt);
          if (!(tmask[jj * LEN + ji] <= 0)) {
              if ((tmask[(jj - 1) * LEN + ji] < 0)) {
                ssha[jj * LEN + ji] = (amp_tide * native_sin((omega_tide * rtime)));
              } else {
                if ((tmask[(jj + 1) * LEN + ji] < 0)) {
                  ssha[jj * LEN + ji] = (amp_tide * native_sin((omega_tide * rtime)));
                } else {
                  if ((tmask[jj * LEN + (ji + 1)] < 0)) {
                    ssha[jj * LEN + ji] = (amp_tide * native_sin((omega_tide * rtime)));
                  } else {
                    if ((tmask[jj * LEN + (ji - 1)] < 0)) {
                      ssha[jj * LEN + ji] = (amp_tide * native_sin((omega_tide * rtime)));
                    }
                  }
                }
              }
          }
      }
  }
}

__attribute__((vec_type_hint(double)))
__kernel __attribute__ ((reqd_work_group_size(1, 1, 1)))
__kernel __attribute__((xcl_zero_global_work_offset))
__kernel void continuity_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict ssha,
  __global double * restrict sshn,
  __global double * restrict sshn_u,
  __global double * restrict sshn_v,
  __global double * restrict hu,
  __global double * restrict hv,
  __global double * restrict un,
  __global double * restrict vn,
  __global double * restrict e12t,
  double rdt
  ){
  for (int jj = ystart; jj <= ystop; jj++){
      __attribute__((opencl_unroll_hint(8)))
      for (int ji = xstart; ji <= xstop; ji++){
          double rtmp1 = ((sshn_u[jj * LEN + ji] + hu[jj * LEN + ji]) * un[jj * LEN + ji]);
          double rtmp2 = ((sshn_u[jj * LEN + (ji - 1)] + hu[jj * LEN + (ji - 1)]) * un[jj * LEN + (ji - 1)]);
          double rtmp3 = ((sshn_v[jj * LEN + ji] + hv[jj * LEN + ji]) * vn[jj * LEN + ji]);
          double rtmp4 = ((sshn_v[(jj - 1) * LEN + ji] + hv[(jj - 1) * LEN + ji]) * vn[(jj - 1) * LEN + ji]);
          ssha[jj * LEN + ji] = (sshn[jj * LEN + ji] + (((((rtmp2 - rtmp1) + rtmp4) - rtmp3) * rdt) / e12t[jj * LEN + ji]));
      }
  }
}

__kernel void field_copy_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict output,
  __global double * restrict input
  ){
  for (int jj = ystart; jj <= ystop; jj++){
      for (int ji = xstart; ji <= xstop; ji++){
        output[jj * LEN + ji] = input[jj * LEN + ji];
      }
  }
}

__kernel void momentum_u_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict ua,
  __global double * restrict un,
  __global double * restrict vn,
  __global double * restrict hu,
  __global double * restrict hv,
  __global double * restrict ht,
  __global double * restrict ssha_u,
  __global double * restrict sshn,
  __global double * restrict sshn_u,
  __global double * restrict sshn_v,
  __global int * restrict tmask,
  __global double * restrict e1u,
  __global double * restrict e1v,
  __global double * restrict e1t,
  __global double * restrict e2u,
  __global double * restrict e2t,
  __global double * restrict e12u,
  __global double * restrict gphiu,
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
  for (int jj = ystart; jj <= ystop; jj++){
      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask[jj * LEN + ji] + tmask[jj * LEN + (ji + 1)]) <= 0)) {
            continue;
          }
          if (((tmask[jj * LEN + ji] <= 0) || (tmask[jj * LEN + (ji + 1)] <= 0))) {
            continue;
          }
          u_e = ((0.5 * (un[jj * LEN + ji] + un[jj * LEN + (ji + 1)])) * e2t[jj * LEN + (ji + 1)]);
          depe = (ht[jj * LEN + (ji + 1)] + sshn[jj * LEN + (ji + 1)]);
          u_w = ((0.5 * (un[jj * LEN + ji] + un[jj * LEN + (ji - 1)])) * e2t[jj * LEN + ji]);
          depw = (ht[jj * LEN + ji] + sshn[jj * LEN + ji]);
          v_sc = (0.5 * (vn[(jj - 1) * LEN + ji] + vn[(jj - 1) * LEN + (ji + 1)]));
          v_s = ((0.5 * v_sc) * (e1v[(jj - 1) * LEN + ji] + e1v[(jj - 1) * LEN + (ji + 1)]));
          deps = (0.5 * (((hv[(jj - 1) * LEN + ji] + sshn_v[(jj - 1) * LEN + ji]) + hv[(jj - 1) * LEN + (ji + 1)]) + sshn_v[(jj - 1) * LEN + (ji + 1)]));
          v_nc = (0.5 * (vn[jj * LEN + ji] + vn[jj * LEN + (ji + 1)]));
          v_n = ((0.5 * v_nc) * (e1v[jj * LEN + ji] + e1v[jj * LEN + (ji + 1)]));
          depn = (0.5 * (((hv[jj * LEN + ji] + sshn_v[jj * LEN + ji]) + hv[jj * LEN + (ji + 1)]) + sshn_v[jj * LEN + (ji + 1)]));
          uu_w = (((0.5 - copysign(0.5, u_w)) * un[jj * LEN + ji]) + ((0.5 + copysign(0.5, u_w)) * un[jj * LEN + (ji - 1)]));
          uu_e = (((0.5 + copysign(0.5, u_e)) * un[jj * LEN + ji]) + ((0.5 - copysign(0.5, u_e)) * un[jj * LEN + (ji + 1)]));
          if (((tmask[(jj - 1) * LEN + ji] <= 0) || (tmask[(jj - 1) * LEN + (ji + 1)] <= 0))) {
            uu_s = ((0.5 - copysign(0.5, v_s)) * un[jj * LEN + ji]);
          } else {
            uu_s = (((0.5 - copysign(0.5, v_s)) * un[jj * LEN + ji]) + ((0.5 + copysign(0.5, v_s)) * un[(jj - 1) * LEN + ji]));
          }
          if (((tmask[(jj + 1) * LEN + ji] <= 0) || (tmask[(jj + 1) * LEN + (ji + 1)] <= 0))) {
            uu_n = ((0.5 + copysign(0.5, v_n)) * un[jj * LEN + ji]);
          } else {
            uu_n = (((0.5 + copysign(0.5, v_n)) * un[jj * LEN + ji]) + ((0.5 - copysign(0.5, v_n)) * un[(jj + 1) * LEN + ji]));
          }
          adv = (((((uu_w * u_w) * depw) - ((uu_e * u_e) * depe)) + ((uu_s * v_s) * deps)) - ((uu_n * v_n) * depn));
          dudx_e = (((un[jj * LEN + (ji + 1)] - un[jj * LEN + ji]) / e1t[jj * LEN + (ji + 1)]) * (ht[jj * LEN + (ji + 1)] + sshn[jj * LEN + (ji + 1)]));
          dudx_w = (((un[jj * LEN + ji] - un[jj * LEN + (ji - 1)]) / e1t[jj * LEN + ji]) * (ht[jj * LEN + ji] + sshn[jj * LEN + ji]));
          if (((tmask[(jj - 1) * LEN + ji] <= 0) || (tmask[(jj - 1) * LEN + (ji + 1)] <= 0))) {
            dudy_s = 0.0;
          } else {
            dudy_s = (((un[jj * LEN + ji] - un[(jj - 1) * LEN + ji]) / (e2u[jj * LEN + ji] + e2u[(jj - 1) * LEN + ji])) * (((hu[jj * LEN + ji] + sshn_u[jj * LEN + ji]) + hu[(jj - 1) * LEN + ji]) + sshn_u[(jj - 1) * LEN + ji]));
          }
          if (((tmask[(jj + 1) * LEN + ji] <= 0) || (tmask[(jj + 1) * LEN + (ji + 1)] <= 0))) {
            dudy_n = 0.0;
          } else {
            dudy_n = (((un[(jj + 1) * LEN + ji] - un[jj * LEN + ji]) / (e2u[jj * LEN + ji] + e2u[(jj + 1) * LEN + ji])) * (((hu[jj * LEN + ji] + sshn_u[jj * LEN + ji]) + hu[(jj + 1) * LEN + ji]) + sshn_u[(jj + 1) * LEN + ji]));
          }
          vis = (((dudx_e - dudx_w) * e2u[jj * LEN + ji]) + (((dudy_n - dudy_s) * e1u[jj * LEN + ji]) * 0.5));
          vis = (visc * vis);
          cor = (((0.5 * (((2. * omega) * sin((gphiu[jj * LEN + ji] * d2r))) * (v_sc + v_nc))) * e12u[jj * LEN + ji]) * (hu[jj * LEN + ji] + sshn_u[jj * LEN + ji]));
          hpg = (-(((g * (hu[jj * LEN + ji] + sshn_u[jj * LEN + ji])) * e2u[jj * LEN + ji]) * (sshn[jj * LEN + (ji + 1)] - sshn[jj * LEN + ji])));
          ua[jj * LEN + ji] = ((((un[jj * LEN + ji] * (hu[jj * LEN + ji] + sshn_u[jj * LEN + ji])) + ((rdt * (((adv + vis) + cor) + hpg)) / e12u[jj * LEN + ji])) / (hu[jj * LEN + ji] + ssha_u[jj * LEN + ji])) / (1.0 + (cbfr * rdt)));
      }
  }
}

__kernel void momentum_v_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict va,
  __global double * restrict un,
  __global double * restrict vn,
  __global double * restrict hu,
  __global double * restrict hv,
  __global double * restrict ht,
  __global double * restrict ssha_v,
  __global double * restrict sshn,
  __global double * restrict sshn_u,
  __global double * restrict sshn_v,
  __global int * restrict tmask,
  __global double * restrict e1v,
  __global double * restrict e1t,
  __global double * restrict e2u,
  __global double * restrict e2v,
  __global double * restrict e2t,
  __global double * restrict e12v,
  __global double * restrict gphiv,
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
  for (int jj = ystart; jj <= ystop; jj++){
      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask[jj * LEN + ji] + tmask[jj * LEN + (ji + 1)]) <= 0)) {
            continue;
          }
          if (((tmask[jj * LEN + ji] <= 0) || (tmask[(jj + 1) * LEN + ji] <= 0))) {
            continue;
          }
          v_n = ((0.5 * (vn[jj * LEN + ji] + vn[(jj + 1) * LEN + ji])) * e1t[(jj + 1) * LEN + ji]);
          depn = (ht[(jj + 1) * LEN + ji] + sshn[(jj + 1) * LEN + ji]);
          v_s = ((0.5 * (vn[jj * LEN + ji] + vn[(jj - 1) * LEN + ji])) * e1t[jj * LEN + ji]);
          deps = (ht[jj * LEN + ji] + sshn[jj * LEN + ji]);
          u_wc = (0.5 * (un[jj * LEN + (ji - 1)] + un[(jj + 1) * LEN + (ji - 1)]));
          u_w = ((0.5 * u_wc) * (e2u[jj * LEN + (ji - 1)] + e2u[(jj + 1) * LEN + (ji - 1)]));
          depw = (0.50 * (((hu[jj * LEN + (ji - 1)] + sshn_u[jj * LEN + (ji - 1)]) + hu[(jj + 1) * LEN + (ji - 1)]) + sshn_u[(jj + 1) * LEN + (ji - 1)]));
          u_ec = (0.5 * (un[jj * LEN + ji] + un[(jj + 1) * LEN + ji]));
          u_e = ((0.5 * u_ec) * (e2u[jj * LEN + ji] + e2u[(jj + 1) * LEN + ji]));
          depe = (0.50 * (((hu[jj * LEN + ji] + sshn_u[jj * LEN + ji]) + hu[(jj + 1) * LEN + ji]) + sshn_u[(jj + 1) * LEN + ji]));
          vv_s = (((0.5 - copysign(0.5, v_s)) * vn[jj * LEN + ji]) + ((0.5 + copysign(0.5, v_s)) * vn[(jj - 1) * LEN + ji]));
          vv_n = (((0.5 + copysign(0.5, v_n)) * vn[jj * LEN + ji]) + ((0.5 - copysign(0.5, v_n)) * vn[(jj + 1) * LEN + ji]));
          if (((tmask[jj * LEN + (ji - 1)] <= 0) || (tmask[(jj + 1) * LEN + (ji - 1)] <= 0))) {
            vv_w = ((0.5 - copysign(0.5, u_w)) * vn[jj * LEN + ji]);
          } else {
            vv_w = (((0.5 - copysign(0.5, u_w)) * vn[jj * LEN + ji]) + ((0.5 + copysign(0.5, u_w)) * vn[jj * LEN + (ji - 1)]));
          }
          if (((tmask[jj * LEN + (ji + 1)] <= 0) || (tmask[(jj + 1) * LEN + (ji + 1)] <= 0))) {
            vv_e = ((0.5 + copysign(0.5, u_e)) * vn[jj * LEN + ji]);
          } else {
            vv_e = (((0.5 + copysign(0.5, u_e)) * vn[jj * LEN + ji]) + ((0.5 - copysign(0.5, u_e)) * vn[jj * LEN + (ji + 1)]));
          }
          adv = (((((vv_w * u_w) * depw) - ((vv_e * u_e) * depe)) + ((vv_s * v_s) * deps)) - ((vv_n * v_n) * depn));
          dvdy_n = (((vn[(jj + 1) * LEN + ji] - vn[jj * LEN + ji]) / e2t[(jj + 1) * LEN + ji]) * (ht[(jj + 1) * LEN + ji] + sshn[(jj + 1) * LEN + ji]));
          dvdy_s = (((vn[jj * LEN + ji] - vn[(jj - 1) * LEN + ji]) / e2t[jj * LEN + ji]) * (ht[jj * LEN + ji] + sshn[jj * LEN + ji]));
          if (((tmask[jj * LEN + (ji - 1)] <= 0) || (tmask[(jj + 1) * LEN + (ji - 1)] <= 0))) {
            dvdx_w = 0.0;
          } else {
            dvdx_w = (((vn[jj * LEN + ji] - vn[jj * LEN + (ji - 1)]) / (e1v[jj * LEN + ji] + e1v[jj * LEN + (ji - 1)])) * (((hv[jj * LEN + ji] + sshn_v[jj * LEN + ji]) + hv[jj * LEN + (ji - 1)]) + sshn_v[jj * LEN + (ji - 1)]));
          }
          if (((tmask[jj * LEN + (ji + 1)] <= 0) || (tmask[(jj + 1) * LEN + (ji + 1)] <= 0))) {
            dvdx_e = 0.0;
          } else {
            dvdx_e = (((vn[jj * LEN + (ji + 1)] - vn[jj * LEN + ji]) / (e1v[jj * LEN + ji] + e1v[jj * LEN + (ji + 1)])) * (((hv[jj * LEN + ji] + sshn_v[jj * LEN + ji]) + hv[jj * LEN + (ji + 1)]) + sshn_v[jj * LEN + (ji + 1)]));
          }
          vis = (((dvdy_n - dvdy_s) * e1v[jj * LEN + ji]) + (((dvdx_e - dvdx_w) * e2v[jj * LEN + ji]) * 0.5));
          vis = (visc * vis);
          cor = (-(((0.5 * (((2. * omega) * sin((gphiv[jj * LEN + ji] * d2r))) * (u_ec + u_wc))) * e12v[jj * LEN + ji]) * (hv[jj * LEN + ji] + sshn_v[jj * LEN + ji])));
          hpg = (-(((g * (hv[jj * LEN + ji] + sshn_v[jj * LEN + ji])) * e1v[jj * LEN + ji]) * (sshn[(jj + 1) * LEN + ji] - sshn[jj * LEN + ji])));
          va[jj * LEN + ji] = ((((vn[jj * LEN + ji] * (hv[jj * LEN + ji] + sshn_v[jj * LEN + ji])) + ((rdt * (((adv + vis) + cor) + hpg)) / e12v[jj * LEN + ji])) / (hv[jj * LEN + ji] + ssha_v[jj * LEN + ji])) / (1.0 + (cbfr * rdt)));
      }
  }
}

__kernel void next_sshu_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict sshn_u,
  __global double * restrict sshn,
  __global int * restrict tmask,
  __global double * restrict e12t,
  __global double * restrict e12u
  ){
  double rtmp1;
  for (int jj = ystart; jj <= ystop; jj++){
      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask[jj * LEN + ji] + tmask[jj * LEN + (ji + 1)]) <= 0)) {
            continue;
          }
          if (((tmask[jj * LEN + ji] * tmask[jj * LEN + (ji + 1)]) > 0)) {
            rtmp1 = ((e12t[jj * LEN + ji] * sshn[jj * LEN + ji]) + (e12t[jj * LEN + (ji + 1)] * sshn[jj * LEN + (ji + 1)]));
            sshn_u[jj * LEN + ji] = ((0.5 * rtmp1) / e12u[jj * LEN + ji]);
          } else {
            if ((tmask[jj * LEN + ji] <= 0)) {
              sshn_u[jj * LEN + ji] = sshn[jj * LEN + (ji + 1)];
            } else {
              if ((tmask[jj * LEN + (ji + 1)] <= 0)) {
                sshn_u[jj * LEN + ji] = sshn[jj * LEN + ji];
              }
            }
          }
      }
  }
}

__kernel void next_sshv_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict sshn_v,
  __global double * restrict sshn,
  __global int * restrict tmask,
  __global double * restrict e12t,
  __global double * restrict e12v
  ){
  double rtmp1;
  for (int jj = ystart; jj <= ystop; jj++){
      for (int ji = xstart; ji <= xstop; ji++){
          if (((tmask[jj * LEN + ji] + tmask[(jj + 1) * LEN + ji]) <= 0)) {
            continue;
          }
          if (((tmask[jj * LEN + ji] * tmask[(jj + 1) * LEN + ji]) > 0)) {
            rtmp1 = ((e12t[jj * LEN + ji] * sshn[jj * LEN + ji]) + (e12t[(jj + 1) * LEN + ji] * sshn[(jj + 1) * LEN + ji]));
            sshn_v[jj * LEN + ji] = ((0.5 * rtmp1) / e12v[jj * LEN + ji]);
          } else {
            if ((tmask[jj * LEN + ji] <= 0)) {
              sshn_v[jj * LEN + ji] = sshn[(jj + 1) * LEN + ji];
            } else {
              if ((tmask[(jj + 1) * LEN + ji] <= 0)) {
                sshn_v[jj * LEN + ji] = sshn[jj * LEN + ji];
              }
            }
          }
      }
  }
}

