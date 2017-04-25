#ifndef __CONTINUITY_HEADER
#define __CONTINUITY_HEADER

void continuity_code(int ji, int jj, int width,                     
		     double *ssha,
		     double *sshn,
		     double *sshn_u,
		     double *sshn_v,
		     double* hu,
		     double *hv,
		     double *un,
		     double *vn,
		     double rdt,
		     double *e12t);
#endif
