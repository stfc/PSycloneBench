#ifndef _TIME_UPDATE_INCLUDE
#define _TIME_UPDATE_INCLUDE

void next_sshu_code(int ji, int jj, int width,
		    double *sshn_u, double *sshn,
		    int *tmask, double *e12t, double *e12u);

void next_sshv_code(int ji, int jj, int width,
		    double *sshn_v, double *sshn, int *tmask,
		    double *e12t, double *e12v);

#endif
