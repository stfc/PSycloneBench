#ifndef _MOM_KERNEL_INCLUDE
#define _MOM_KERNEL_INCLUDE
void momentum_u_code(int ji, int jj, int width,
		     double *ua, double *un, double *vn,
		     double *hu, double *hv, double *ht, double *ssha_u,
		     double *sshn, double *sshn_u, double *sshn_v,
		     int *tmask,
		     double *e1u, double *e1v, double *e1t,
		     double *e2u, double *e2t, double *e12u, double *gphiu,
		     double rdt, double cbfr, double visc);

void momentum_v_code(int ji, int jj, int width,
		     double *va, double *un, double *vn, 
		     double *hu, double *hv, double *ht, double *ssha_v, 
		     double *sshn, double *sshn_u, double *sshn_v, 
		     int *tmask, double *e1v, double *e1t, double *e2u,
		     double *e2v, double *e2t, double *e12v, double *gphiv,
		     double rdt, double cbfr, double visc);
#endif
