#include <iostream>
#include <chrono>
#include <vector>
#include <cstdlib>

// Compared to the Fortran implementations:
//  - Swap loop order

extern "C" void c_invoke_time_step(
        // Fields
        double * ssha_t,
        double * sshn_t,
        double * sshn_u,
        double * sshn_v,
        double * hu,
        double * hv,
        double * un,
        double * vn,
        double * ua,
        double * ht,
        double * ssha_u,
        double * va,
        double * ssha_v,
        // Grid
        double * e12t,
        // Scalars
        int istp,
        int nx,
        int ny,
        double rdt
        ){

    int M = nx;
    int N = ny;
    for(int jj = 20; jj < M - 20; jj++){
        for(int ji = 20; ji < N - 20; ji++){
            // Continuity kernel
            int idx = jj * M + ji;
            int idxim1 = idx - 1;
            int idxjm1 = idx - M;

            double rtmp1 = (sshn_u[idx]    + hu[idx])    * un[idx];
            double rtmp2 = (sshn_u[idxim1] + hu[idxim1]) * un[idxim1];
            double rtmp3 = (sshn_v[idx]    + hv[idx])    * vn[idx];
            double rtmp4 = (sshn_v[idxjm1] + hv[idxjm1]) * vn[idxjm1];

            ssha_t[idx] = sshn_t[idx] + (rtmp2 - rtmp1 + rtmp4 - rtmp3) *
                rdt / e12t[idx];

        }
    }

    
}
