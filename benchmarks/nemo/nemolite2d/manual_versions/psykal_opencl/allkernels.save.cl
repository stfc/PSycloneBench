__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
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
  int jiu;
  int uaLEN1 = get_global_size(0);
  int uaLEN2 = get_global_size(1);
  int huLEN1 = get_global_size(0);
  int huLEN2 = get_global_size(1);
  int sshn_uLEN1 = get_global_size(0);
  int sshn_uLEN2 = get_global_size(1);
  int tmaskLEN1 = get_global_size(0);
  int tmaskLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] + tmask[jj * tmaskLEN1 + (ji + 1)]) <= (-1))) {
    return;
  }
  if ((tmask[jj * tmaskLEN1 + ji] < 0)) {
    jiu = (ji + 1);
    ua[jj * uaLEN1 + ji] = (ua[jj * uaLEN1 + jiu] + (native_sqrt((g / hu[jj * huLEN1 + ji])) * (sshn_u[jj * sshn_uLEN1 + ji] - sshn_u[jj * sshn_uLEN1 + jiu])));
  } else {
    if ((tmask[jj * tmaskLEN1 + (ji + 1)] < 0)) {
      jiu = (ji - 1);
      ua[jj * uaLEN1 + ji] = (ua[jj * uaLEN1 + jiu] + (native_sqrt((g / hu[jj * huLEN1 + ji])) * (sshn_u[jj * sshn_uLEN1 + ji] - sshn_u[jj * sshn_uLEN1 + jiu])));
    }
  }
}

__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
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
  int jiv;
  int vaLEN1 = get_global_size(0);
  int vaLEN2 = get_global_size(1);
  int hvLEN1 = get_global_size(0);
  int hvLEN2 = get_global_size(1);
  int sshn_vLEN1 = get_global_size(0);
  int sshn_vLEN2 = get_global_size(1);
  int tmaskLEN1 = get_global_size(0);
  int tmaskLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] + tmask[(jj + 1) * tmaskLEN1 + ji]) <= (-1))) {
    return;
  }
  if ((tmask[jj * tmaskLEN1 + ji] < 0)) {
    jiv = (jj + 1);
    va[jj * vaLEN1 + ji] = (va[jiv * vaLEN1 + ji] + (native_sqrt((g / hv[jj * hvLEN1 + ji])) * (sshn_v[jj * sshn_vLEN1 + ji] - sshn_v[jiv * sshn_vLEN1 + ji])));
  } else {
    if ((tmask[(jj + 1) * tmaskLEN1 + ji] < 0)) {
      jiv = (jj - 1);
      va[jj * vaLEN1 + ji] = (va[jiv * vaLEN1 + ji] + (native_sqrt((g / hv[jj * hvLEN1 + ji])) * (sshn_v[jj * sshn_vLEN1 + ji] - sshn_v[jiv * sshn_vLEN1 + ji])));
    }
  }
}

__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
__attribute__((xcl_zero_global_work_offset))
__kernel void bc_solid_u_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict ua,
  const __global int * restrict tmask
  ){
  int uaLEN1 = get_global_size(0);
  int uaLEN2 = get_global_size(1);
  int tmaskLEN1 = get_global_size(0);
  int tmaskLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] * tmask[jj * tmaskLEN1 + (ji + 1)]) == 0)) {
    ua[jj * uaLEN1 + ji] = 0.;
  }
}

__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
__attribute__((xcl_zero_global_work_offset))
__kernel void bc_solid_v_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict va,
  const __global int * restrict tmask
  ){
  int vaLEN1 = get_global_size(0);
  int vaLEN2 = get_global_size(1);
  int tmaskLEN1 = get_global_size(0);
  int tmaskLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] * tmask[(jj + 1) * tmaskLEN1 + ji]) == 0)) {
    va[jj * vaLEN1 + ji] = 0.;
  }
}

__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
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
  double amp_tide;
  double omega_tide;
  double rtime;
  int sshaLEN1 = get_global_size(0);
  int sshaLEN2 = get_global_size(1);
  int tmaskLEN1 = get_global_size(0);
  int tmaskLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  amp_tide = 0.2;
  omega_tide = ((2.0 * 3.14159) / (12.42 * 3600.));
  rtime = ((float)istep * rdt);
  if ((tmask[jj * tmaskLEN1 + ji] <= 0)) {
    return;
  }
  if ((tmask[(jj - 1) * tmaskLEN1 + ji] < 0)) {
    ssha[jj * sshaLEN1 + ji] = (amp_tide * native_sin((omega_tide * rtime)));
  } else {
    if ((tmask[(jj + 1) * tmaskLEN1 + ji] < 0)) {
      ssha[jj * sshaLEN1 + ji] = (amp_tide * native_sin((omega_tide * rtime)));
    } else {
      if ((tmask[jj * tmaskLEN1 + (ji + 1)] < 0)) {
        ssha[jj * sshaLEN1 + ji] = (amp_tide * native_sin((omega_tide * rtime)));
      } else {
        if ((tmask[jj * tmaskLEN1 + (ji - 1)] < 0)) {
          ssha[jj * sshaLEN1 + ji] = (amp_tide * native_sin((omega_tide * rtime)));
        }
      }
    }
  }
}

__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
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
  double rtmp1;
  double rtmp2;
  double rtmp3;
  double rtmp4;
  int sshaLEN1 = get_global_size(0);
  int sshaLEN2 = get_global_size(1);
  int sshnLEN1 = get_global_size(0);
  int sshnLEN2 = get_global_size(1);
  int sshn_uLEN1 = get_global_size(0);
  int sshn_uLEN2 = get_global_size(1);
  int sshn_vLEN1 = get_global_size(0);
  int sshn_vLEN2 = get_global_size(1);
  int huLEN1 = get_global_size(0);
  int huLEN2 = get_global_size(1);
  int hvLEN1 = get_global_size(0);
  int hvLEN2 = get_global_size(1);
  int unLEN1 = get_global_size(0);
  int unLEN2 = get_global_size(1);
  int vnLEN1 = get_global_size(0);
  int vnLEN2 = get_global_size(1);
  int e12tLEN1 = get_global_size(0);
  int e12tLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  rtmp1 = ((sshn_u[jj * sshn_uLEN1 + ji] + hu[jj * huLEN1 + ji]) * un[jj * unLEN1 + ji]);
  rtmp2 = ((sshn_u[jj * sshn_uLEN1 + (ji - 1)] + hu[jj * huLEN1 + (ji - 1)]) * un[jj * unLEN1 + (ji - 1)]);
  rtmp3 = ((sshn_v[jj * sshn_vLEN1 + ji] + hv[jj * hvLEN1 + ji]) * vn[jj * vnLEN1 + ji]);
  rtmp4 = ((sshn_v[(jj - 1) * sshn_vLEN1 + ji] + hv[(jj - 1) * hvLEN1 + ji]) * vn[(jj - 1) * vnLEN1 + ji]);
  ssha[jj * sshaLEN1 + ji] = (sshn[jj * sshnLEN1 + ji] + (((((rtmp2 - rtmp1) + rtmp4) - rtmp3) * rdt) / e12t[jj * e12tLEN1 + ji]));
}

__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
__attribute__((xcl_zero_global_work_offset))
__kernel void field_copy_code(
  int xstart,
  int xstop,
  int ystart,
  int ystop,
  __global double * restrict output,
  const __global double * restrict input
  ){
  int outputLEN1 = get_global_size(0);
  int outputLEN2 = get_global_size(1);
  int inputLEN1 = get_global_size(0);
  int inputLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  output[jj * outputLEN1 + ji] = input[jj * inputLEN1 + ji];
}

__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
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
  int uaLEN1 = get_global_size(0);
  int uaLEN2 = get_global_size(1);
  int unLEN1 = get_global_size(0);
  int unLEN2 = get_global_size(1);
  int vnLEN1 = get_global_size(0);
  int vnLEN2 = get_global_size(1);
  int huLEN1 = get_global_size(0);
  int huLEN2 = get_global_size(1);
  int hvLEN1 = get_global_size(0);
  int hvLEN2 = get_global_size(1);
  int htLEN1 = get_global_size(0);
  int htLEN2 = get_global_size(1);
  int ssha_uLEN1 = get_global_size(0);
  int ssha_uLEN2 = get_global_size(1);
  int sshnLEN1 = get_global_size(0);
  int sshnLEN2 = get_global_size(1);
  int sshn_uLEN1 = get_global_size(0);
  int sshn_uLEN2 = get_global_size(1);
  int sshn_vLEN1 = get_global_size(0);
  int sshn_vLEN2 = get_global_size(1);
  int tmaskLEN1 = get_global_size(0);
  int tmaskLEN2 = get_global_size(1);
  int e1uLEN1 = get_global_size(0);
  int e1uLEN2 = get_global_size(1);
  int e1vLEN1 = get_global_size(0);
  int e1vLEN2 = get_global_size(1);
  int e1tLEN1 = get_global_size(0);
  int e1tLEN2 = get_global_size(1);
  int e2uLEN1 = get_global_size(0);
  int e2uLEN2 = get_global_size(1);
  int e2tLEN1 = get_global_size(0);
  int e2tLEN2 = get_global_size(1);
  int e12uLEN1 = get_global_size(0);
  int e12uLEN2 = get_global_size(1);
  int gphiuLEN1 = get_global_size(0);
  int gphiuLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] + tmask[jj * tmaskLEN1 + (ji + 1)]) <= 0)) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] <= 0) || (tmask[jj * tmaskLEN1 + (ji + 1)] <= 0))) {
    return;
  }
  u_e = ((0.5 * (un[jj * unLEN1 + ji] + un[jj * unLEN1 + (ji + 1)])) * e2t[jj * e2tLEN1 + (ji + 1)]);
  depe = (ht[jj * htLEN1 + (ji + 1)] + sshn[jj * sshnLEN1 + (ji + 1)]);
  u_w = ((0.5 * (un[jj * unLEN1 + ji] + un[jj * unLEN1 + (ji - 1)])) * e2t[jj * e2tLEN1 + ji]);
  depw = (ht[jj * htLEN1 + ji] + sshn[jj * sshnLEN1 + ji]);
  v_sc = (0.5 * (vn[(jj - 1) * vnLEN1 + ji] + vn[(jj - 1) * vnLEN1 + (ji + 1)]));
  v_s = ((0.5 * v_sc) * (e1v[(jj - 1) * e1vLEN1 + ji] + e1v[(jj - 1) * e1vLEN1 + (ji + 1)]));
  deps = (0.5 * (((hv[(jj - 1) * hvLEN1 + ji] + sshn_v[(jj - 1) * sshn_vLEN1 + ji]) + hv[(jj - 1) * hvLEN1 + (ji + 1)]) + sshn_v[(jj - 1) * sshn_vLEN1 + (ji + 1)]));
  v_nc = (0.5 * (vn[jj * vnLEN1 + ji] + vn[jj * vnLEN1 + (ji + 1)]));
  v_n = ((0.5 * v_nc) * (e1v[jj * e1vLEN1 + ji] + e1v[jj * e1vLEN1 + (ji + 1)]));
  depn = (0.5 * (((hv[jj * hvLEN1 + ji] + sshn_v[jj * sshn_vLEN1 + ji]) + hv[jj * hvLEN1 + (ji + 1)]) + sshn_v[jj * sshn_vLEN1 + (ji + 1)]));
  uu_w = (((0.5 - copysign(0.5, u_w)) * un[jj * unLEN1 + ji]) + ((0.5 + copysign(0.5, u_w)) * un[jj * unLEN1 + (ji - 1)]));
  uu_e = (((0.5 + copysign(0.5, u_e)) * un[jj * unLEN1 + ji]) + ((0.5 - copysign(0.5, u_e)) * un[jj * unLEN1 + (ji + 1)]));
  if (((tmask[(jj - 1) * tmaskLEN1 + ji] <= 0) || (tmask[(jj - 1) * tmaskLEN1 + (ji + 1)] <= 0))) {
    uu_s = ((0.5 - copysign(0.5, v_s)) * un[jj * unLEN1 + ji]);
  } else {
    uu_s = (((0.5 - copysign(0.5, v_s)) * un[jj * unLEN1 + ji]) + ((0.5 + copysign(0.5, v_s)) * un[(jj - 1) * unLEN1 + ji]));
  }
  if (((tmask[(jj + 1) * tmaskLEN1 + ji] <= 0) || (tmask[(jj + 1) * tmaskLEN1 + (ji + 1)] <= 0))) {
    uu_n = ((0.5 + copysign(0.5, v_n)) * un[jj * unLEN1 + ji]);
  } else {
    uu_n = (((0.5 + copysign(0.5, v_n)) * un[jj * unLEN1 + ji]) + ((0.5 - copysign(0.5, v_n)) * un[(jj + 1) * unLEN1 + ji]));
  }
  adv = (((((uu_w * u_w) * depw) - ((uu_e * u_e) * depe)) + ((uu_s * v_s) * deps)) - ((uu_n * v_n) * depn));
  dudx_e = (((un[jj * unLEN1 + (ji + 1)] - un[jj * unLEN1 + ji]) / e1t[jj * e1tLEN1 + (ji + 1)]) * (ht[jj * htLEN1 + (ji + 1)] + sshn[jj * sshnLEN1 + (ji + 1)]));
  dudx_w = (((un[jj * unLEN1 + ji] - un[jj * unLEN1 + (ji - 1)]) / e1t[jj * e1tLEN1 + ji]) * (ht[jj * htLEN1 + ji] + sshn[jj * sshnLEN1 + ji]));
  if (((tmask[(jj - 1) * tmaskLEN1 + ji] <= 0) || (tmask[(jj - 1) * tmaskLEN1 + (ji + 1)] <= 0))) {
    dudy_s = 0.0;
  } else {
    dudy_s = (((un[jj * unLEN1 + ji] - un[(jj - 1) * unLEN1 + ji]) / (e2u[jj * e2uLEN1 + ji] + e2u[(jj - 1) * e2uLEN1 + ji])) * (((hu[jj * huLEN1 + ji] + sshn_u[jj * sshn_uLEN1 + ji]) + hu[(jj - 1) * huLEN1 + ji]) + sshn_u[(jj - 1) * sshn_uLEN1 + ji]));
  }
  if (((tmask[(jj + 1) * tmaskLEN1 + ji] <= 0) || (tmask[(jj + 1) * tmaskLEN1 + (ji + 1)] <= 0))) {
    dudy_n = 0.0;
  } else {
    dudy_n = (((un[(jj + 1) * unLEN1 + ji] - un[jj * unLEN1 + ji]) / (e2u[jj * e2uLEN1 + ji] + e2u[(jj + 1) * e2uLEN1 + ji])) * (((hu[jj * huLEN1 + ji] + sshn_u[jj * sshn_uLEN1 + ji]) + hu[(jj + 1) * huLEN1 + ji]) + sshn_u[(jj + 1) * sshn_uLEN1 + ji]));
  }
  vis = (((dudx_e - dudx_w) * e2u[jj * e2uLEN1 + ji]) + (((dudy_n - dudy_s) * e1u[jj * e1uLEN1 + ji]) * 0.5));
  vis = (visc * vis);
  cor = (((0.5 * (((2. * omega) * native_sin((gphiu[jj * gphiuLEN1 + ji] * d2r))) * (v_sc + v_nc))) * e12u[jj * e12uLEN1 + ji]) * (hu[jj * huLEN1 + ji] + sshn_u[jj * sshn_uLEN1 + ji]));
  hpg = (-(((g * (hu[jj * huLEN1 + ji] + sshn_u[jj * sshn_uLEN1 + ji])) * e2u[jj * e2uLEN1 + ji]) * (sshn[jj * sshnLEN1 + (ji + 1)] - sshn[jj * sshnLEN1 + ji])));
  ua[jj * uaLEN1 + ji] = ((((un[jj * unLEN1 + ji] * (hu[jj * huLEN1 + ji] + sshn_u[jj * sshn_uLEN1 + ji])) + ((rdt * (((adv + vis) + cor) + hpg)) / e12u[jj * e12uLEN1 + ji])) / (hu[jj * huLEN1 + ji] + ssha_u[jj * ssha_uLEN1 + ji])) / (1.0 + (cbfr * rdt)));
}

__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
__attribute__((xcl_zero_global_work_offset))
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
  int vaLEN1 = get_global_size(0);
  int vaLEN2 = get_global_size(1);
  int unLEN1 = get_global_size(0);
  int unLEN2 = get_global_size(1);
  int vnLEN1 = get_global_size(0);
  int vnLEN2 = get_global_size(1);
  int huLEN1 = get_global_size(0);
  int huLEN2 = get_global_size(1);
  int hvLEN1 = get_global_size(0);
  int hvLEN2 = get_global_size(1);
  int htLEN1 = get_global_size(0);
  int htLEN2 = get_global_size(1);
  int ssha_vLEN1 = get_global_size(0);
  int ssha_vLEN2 = get_global_size(1);
  int sshnLEN1 = get_global_size(0);
  int sshnLEN2 = get_global_size(1);
  int sshn_uLEN1 = get_global_size(0);
  int sshn_uLEN2 = get_global_size(1);
  int sshn_vLEN1 = get_global_size(0);
  int sshn_vLEN2 = get_global_size(1);
  int tmaskLEN1 = get_global_size(0);
  int tmaskLEN2 = get_global_size(1);
  int e1vLEN1 = get_global_size(0);
  int e1vLEN2 = get_global_size(1);
  int e1tLEN1 = get_global_size(0);
  int e1tLEN2 = get_global_size(1);
  int e2uLEN1 = get_global_size(0);
  int e2uLEN2 = get_global_size(1);
  int e2vLEN1 = get_global_size(0);
  int e2vLEN2 = get_global_size(1);
  int e2tLEN1 = get_global_size(0);
  int e2tLEN2 = get_global_size(1);
  int e12vLEN1 = get_global_size(0);
  int e12vLEN2 = get_global_size(1);
  int gphivLEN1 = get_global_size(0);
  int gphivLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] + tmask[jj * tmaskLEN1 + (ji + 1)]) <= 0)) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] <= 0) || (tmask[(jj + 1) * tmaskLEN1 + ji] <= 0))) {
    return;
  }
  v_n = ((0.5 * (vn[jj * vnLEN1 + ji] + vn[(jj + 1) * vnLEN1 + ji])) * e1t[(jj + 1) * e1tLEN1 + ji]);
  depn = (ht[(jj + 1) * htLEN1 + ji] + sshn[(jj + 1) * sshnLEN1 + ji]);
  v_s = ((0.5 * (vn[jj * vnLEN1 + ji] + vn[(jj - 1) * vnLEN1 + ji])) * e1t[jj * e1tLEN1 + ji]);
  deps = (ht[jj * htLEN1 + ji] + sshn[jj * sshnLEN1 + ji]);
  u_wc = (0.5 * (un[jj * unLEN1 + (ji - 1)] + un[(jj + 1) * unLEN1 + (ji - 1)]));
  u_w = ((0.5 * u_wc) * (e2u[jj * e2uLEN1 + (ji - 1)] + e2u[(jj + 1) * e2uLEN1 + (ji - 1)]));
  depw = (0.50 * (((hu[jj * huLEN1 + (ji - 1)] + sshn_u[jj * sshn_uLEN1 + (ji - 1)]) + hu[(jj + 1) * huLEN1 + (ji - 1)]) + sshn_u[(jj + 1) * sshn_uLEN1 + (ji - 1)]));
  u_ec = (0.5 * (un[jj * unLEN1 + ji] + un[(jj + 1) * unLEN1 + ji]));
  u_e = ((0.5 * u_ec) * (e2u[jj * e2uLEN1 + ji] + e2u[(jj + 1) * e2uLEN1 + ji]));
  depe = (0.50 * (((hu[jj * huLEN1 + ji] + sshn_u[jj * sshn_uLEN1 + ji]) + hu[(jj + 1) * huLEN1 + ji]) + sshn_u[(jj + 1) * sshn_uLEN1 + ji]));
  vv_s = (((0.5 - copysign(0.5, v_s)) * vn[jj * vnLEN1 + ji]) + ((0.5 + copysign(0.5, v_s)) * vn[(jj - 1) * vnLEN1 + ji]));
  vv_n = (((0.5 + copysign(0.5, v_n)) * vn[jj * vnLEN1 + ji]) + ((0.5 - copysign(0.5, v_n)) * vn[(jj + 1) * vnLEN1 + ji]));
  if (((tmask[jj * tmaskLEN1 + (ji - 1)] <= 0) || (tmask[(jj + 1) * tmaskLEN1 + (ji - 1)] <= 0))) {
    vv_w = ((0.5 - copysign(0.5, u_w)) * vn[jj * vnLEN1 + ji]);
  } else {
    vv_w = (((0.5 - copysign(0.5, u_w)) * vn[jj * vnLEN1 + ji]) + ((0.5 + copysign(0.5, u_w)) * vn[jj * vnLEN1 + (ji - 1)]));
  }
  if (((tmask[jj * tmaskLEN1 + (ji + 1)] <= 0) || (tmask[(jj + 1) * tmaskLEN1 + (ji + 1)] <= 0))) {
    vv_e = ((0.5 + copysign(0.5, u_e)) * vn[jj * vnLEN1 + ji]);
  } else {
    vv_e = (((0.5 + copysign(0.5, u_e)) * vn[jj * vnLEN1 + ji]) + ((0.5 - copysign(0.5, u_e)) * vn[jj * vnLEN1 + (ji + 1)]));
  }
  adv = (((((vv_w * u_w) * depw) - ((vv_e * u_e) * depe)) + ((vv_s * v_s) * deps)) - ((vv_n * v_n) * depn));
  dvdy_n = (((vn[(jj + 1) * vnLEN1 + ji] - vn[jj * vnLEN1 + ji]) / e2t[(jj + 1) * e2tLEN1 + ji]) * (ht[(jj + 1) * htLEN1 + ji] + sshn[(jj + 1) * sshnLEN1 + ji]));
  dvdy_s = (((vn[jj * vnLEN1 + ji] - vn[(jj - 1) * vnLEN1 + ji]) / e2t[jj * e2tLEN1 + ji]) * (ht[jj * htLEN1 + ji] + sshn[jj * sshnLEN1 + ji]));
  if (((tmask[jj * tmaskLEN1 + (ji - 1)] <= 0) || (tmask[(jj + 1) * tmaskLEN1 + (ji - 1)] <= 0))) {
    dvdx_w = 0.0;
  } else {
    dvdx_w = (((vn[jj * vnLEN1 + ji] - vn[jj * vnLEN1 + (ji - 1)]) / (e1v[jj * e1vLEN1 + ji] + e1v[jj * e1vLEN1 + (ji - 1)])) * (((hv[jj * hvLEN1 + ji] + sshn_v[jj * sshn_vLEN1 + ji]) + hv[jj * hvLEN1 + (ji - 1)]) + sshn_v[jj * sshn_vLEN1 + (ji - 1)]));
  }
  if (((tmask[jj * tmaskLEN1 + (ji + 1)] <= 0) || (tmask[(jj + 1) * tmaskLEN1 + (ji + 1)] <= 0))) {
    dvdx_e = 0.0;
  } else {
    dvdx_e = (((vn[jj * vnLEN1 + (ji + 1)] - vn[jj * vnLEN1 + ji]) / (e1v[jj * e1vLEN1 + ji] + e1v[jj * e1vLEN1 + (ji + 1)])) * (((hv[jj * hvLEN1 + ji] + sshn_v[jj * sshn_vLEN1 + ji]) + hv[jj * hvLEN1 + (ji + 1)]) + sshn_v[jj * sshn_vLEN1 + (ji + 1)]));
  }
  vis = (((dvdy_n - dvdy_s) * e1v[jj * e1vLEN1 + ji]) + (((dvdx_e - dvdx_w) * e2v[jj * e2vLEN1 + ji]) * 0.5));
  vis = (visc * vis);
  cor = (-(((0.5 * (((2. * omega) * native_sin((gphiv[jj * gphivLEN1 + ji] * d2r))) * (u_ec + u_wc))) * e12v[jj * e12vLEN1 + ji]) * (hv[jj * hvLEN1 + ji] + sshn_v[jj * sshn_vLEN1 + ji])));
  hpg = (-(((g * (hv[jj * hvLEN1 + ji] + sshn_v[jj * sshn_vLEN1 + ji])) * e1v[jj * e1vLEN1 + ji]) * (sshn[(jj + 1) * sshnLEN1 + ji] - sshn[jj * sshnLEN1 + ji])));
  va[jj * vaLEN1 + ji] = ((((vn[jj * vnLEN1 + ji] * (hv[jj * hvLEN1 + ji] + sshn_v[jj * sshn_vLEN1 + ji])) + ((rdt * (((adv + vis) + cor) + hpg)) / e12v[jj * e12vLEN1 + ji])) / (hv[jj * hvLEN1 + ji] + ssha_v[jj * ssha_vLEN1 + ji])) / (1.0 + (cbfr * rdt)));
}

__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
__attribute__((xcl_zero_global_work_offset))
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
  int sshn_uLEN1 = get_global_size(0);
  int sshn_uLEN2 = get_global_size(1);
  int sshnLEN1 = get_global_size(0);
  int sshnLEN2 = get_global_size(1);
  int tmaskLEN1 = get_global_size(0);
  int tmaskLEN2 = get_global_size(1);
  int e12tLEN1 = get_global_size(0);
  int e12tLEN2 = get_global_size(1);
  int e12uLEN1 = get_global_size(0);
  int e12uLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] + tmask[jj * tmaskLEN1 + (ji + 1)]) <= 0)) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] * tmask[jj * tmaskLEN1 + (ji + 1)]) > 0)) {
    rtmp1 = ((e12t[jj * e12tLEN1 + ji] * sshn[jj * sshnLEN1 + ji]) + (e12t[jj * e12tLEN1 + (ji + 1)] * sshn[jj * sshnLEN1 + (ji + 1)]));
    sshn_u[jj * sshn_uLEN1 + ji] = ((0.5 * rtmp1) / e12u[jj * e12uLEN1 + ji]);
  } else {
    if ((tmask[jj * tmaskLEN1 + ji] <= 0)) {
      sshn_u[jj * sshn_uLEN1 + ji] = sshn[jj * sshnLEN1 + (ji + 1)];
    } else {
      if ((tmask[jj * tmaskLEN1 + (ji + 1)] <= 0)) {
        sshn_u[jj * sshn_uLEN1 + ji] = sshn[jj * sshnLEN1 + ji];
      }
    }
  }
}

__attribute__((reqd_work_group_size(64, 1, 1)))
__attribute__((vec_type_hint(double)))
__attribute__((xcl_zero_global_work_offset))
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
  int sshn_vLEN1 = get_global_size(0);
  int sshn_vLEN2 = get_global_size(1);
  int sshnLEN1 = get_global_size(0);
  int sshnLEN2 = get_global_size(1);
  int tmaskLEN1 = get_global_size(0);
  int tmaskLEN2 = get_global_size(1);
  int e12tLEN1 = get_global_size(0);
  int e12tLEN2 = get_global_size(1);
  int e12vLEN1 = get_global_size(0);
  int e12vLEN2 = get_global_size(1);
  int ji = get_global_id(0);
  int jj = get_global_id(1);
  if ((((ji < xstart) || (ji > xstop)) || ((jj < ystart) || (jj > ystop)))) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] + tmask[(jj + 1) * tmaskLEN1 + ji]) <= 0)) {
    return;
  }
  if (((tmask[jj * tmaskLEN1 + ji] * tmask[(jj + 1) * tmaskLEN1 + ji]) > 0)) {
    rtmp1 = ((e12t[jj * e12tLEN1 + ji] * sshn[jj * sshnLEN1 + ji]) + (e12t[(jj + 1) * e12tLEN1 + ji] * sshn[(jj + 1) * sshnLEN1 + ji]));
    sshn_v[jj * sshn_vLEN1 + ji] = ((0.5 * rtmp1) / e12v[jj * e12vLEN1 + ji]);
  } else {
    if ((tmask[jj * tmaskLEN1 + ji] <= 0)) {
      sshn_v[jj * sshn_vLEN1 + ji] = sshn[(jj + 1) * sshnLEN1 + ji];
    } else {
      if ((tmask[(jj + 1) * tmaskLEN1 + ji] <= 0)) {
        sshn_v[jj * sshn_vLEN1 + ji] = sshn[jj * sshnLEN1 + ji];
      }
    }
  }
}

