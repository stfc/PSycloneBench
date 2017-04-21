#ifndef _BOUNDARY_CONDITIONS_INCLUDE
#define _BOUNDARY_CONDITIONS_INCLUDE

void bc_ssh_code(int ji, int jj, int width,
		 int istep, double *ssha, int *tmask, double rdt);

void bc_solid_u_code(int ji, int jj, int width, double *ua, int *tmask);
void bc_solid_v_code(int ji, int jj, int width, double *va, int *tmask);
void bc_flather_u_code(int ji, int jj, int width,
		       double *ua, double *hu, double *sshn_u, int *tmask);
void bc_flather_v_code(int ji, int jj, int width,
		       double *va, double *hv, double *sshn_v, int *tmask);

#endif
