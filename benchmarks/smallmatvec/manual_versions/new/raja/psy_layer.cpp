#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "RAJA/RAJA.hpp"
#include "chai/ManagedArray.hpp"


RAJA_DEVICE void matrix_vector_code_kouter(
        int cell, int nlayers, double *lhs, double *x, int ncell_3d, double *matrix, int ndf1, int undf1,
        int *map1, int ndf2, int undf2, int *map2
){

    int ik = cell*nlayers;
    for (int df2 = 0; df2 < ndf2; df2 ++){
        int m2 = map2[df2] - 1; // -1 because map2 contains fortran 1-indexing references
        for (int df = 0; df < ndf1; df ++){
            int m1 = map1[df] - 1; // -1 because map2 contains fortran 1-indexing references
            for (int k = 0; k < nlayers; k ++){
                lhs[m1+k] = lhs[m1+k] + matrix[(ik+k)*ndf1*ndf2 + df2*ndf1 + df] * x[m2+k];
            }
        }
    }
}

RAJA_DEVICE void matrix_vector_code_kouter_atomics(
        int cell, int nlayers, double *lhs, double *x, int ncell_3d, double *matrix, int ndf1, int undf1,
        int *map1, int ndf2, int undf2, int *map2
){

    int ik = cell*nlayers;
    for (int df2 = 0; df2 < ndf2; df2 ++){
        int m2 = map2[df2] - 1; // -1 because map2 contains fortran 1-indexing references
        for (int df = 0; df < ndf1; df ++){
            int m1 = map1[df] - 1; // -1 because map2 contains fortran 1-indexing references
            for (int k = 0; k < nlayers; k ++){    
		RAJA::atomicAdd<RAJA::omp_atomic>(&lhs[m1+k], matrix[(ik+k)*ndf1*ndf2 + df2*ndf1 + df] * x[m2+k]);
            }
        }
    }
}

RAJA_DEVICE void matrix_vector_code_kinner(
        int cell, int nlayers, double *lhs, double *x, int ncell_3d,
        double *matrix, int ndf1, int undf1, int *map1, int ndf2, int undf2, int *map2
){

    for (int df2 = 0; df2 < ndf2; df2 ++){
        int m2 = map2[df2] - 1; // -1 because map2 contains fortran 1-indexing references
        for (int df = 0; df < ndf1; df ++){
            int m1 = map1[df] - 1; // -1 because map2 contains fortran 1-indexing references
            for (int k = 0; k < nlayers; k++){
                 lhs[m1+k]= lhs[m1+k] + matrix[k + cell*nlayers*ndf1*ndf2 + df2*nlayers*ndf1 + df*nlayers] * x[m2+k]; 
            }
        }
    }
}

RAJA_DEVICE void matrix_vector_code_kinner_atomics(
        int cell, int nlayers, double *lhs, double *x, int ncell_3d,
        double *matrix, int ndf1, int undf1, int *map1, int ndf2, int undf2, int *map2
        ){

    for (int df2 = 0; df2 < ndf2; df2 ++){
        int m2 = map2[df2] - 1; // -1 because map2 contains fortran 1-indexing references
        for (int df = 0; df < ndf1; df ++){
            int m1 = map1[df] - 1; // -1 because map2 contains fortran 1-indexing references
            for (int k = 0; k < nlayers; k++){
		RAJA::atomicAdd<RAJA::omp_atomic>(&lhs[m1+k], matrix[k + cell*nlayers*ndf1*ndf2 + df2*nlayers*ndf1 + df*nlayers] * x[m2+k]);
            }
        }
    }
}

extern "C" void c_psy_layer(
        char *traverse, int niters, int ncell, int nlayers,
        int ncell_3d, double *lhs, int *map_lhs,
        int ndf_lhs, int undf_lhs, double *matrix, double *matrix_kinner,
        double *x, int *map_x, int ndf_x, int undf_x,
        int ncolour, int *ncp_colour, int *cmap
) {

    printf("CPP Version\n");

#if defined(USE_CUDA_POLICY)
    using policy = RAJA::cuda_exec<256>;
    printf("Using RAJA cuda version \n");
#elif defined(USE_OPENMP_POLICY)
    using policy = RAJA::omp_parallel_for_exec;
    printf("Using RAJA OpenMP version \n");
#else
    using policy = RAJA::loop_exec;
    printf("Using RAJA sequential version \n");
#endif

    if (memcmp(traverse,"linear-kinner",11)==0){
        printf("Starting computation with linear and kinner\n");
        for (int iter = 1; iter <= niters; iter ++){
            // for (int cell=0; cell < ncell; cell ++){
            //     matrix_vector_code_kinner_atomics(cell, nlayers, lhs, x, ncell_3d, matrix_kinner,
            //             ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
            // }
            RAJA::forall<policy>(RAJA::RangeSegment(0, ncell), [=] RAJA_DEVICE (int cell) {
                 matrix_vector_code_kinner_atomics(cell, nlayers, lhs, x, ncell_3d, matrix_kinner,
                     ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
            });
        }
    }
    else if (memcmp(traverse,"linear",5)==0){
        printf("Linear traversing Version\n");
        for (int iter=1; iter <= niters; iter ++){
            // for (int cell=0; cell < ncell; cell ++){
            //     matrix_vector_code_kouter_atomics(cell, nlayers, lhs, x, ncell_3d, matrix,
            //             ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
            // }
            RAJA::forall<policy>(RAJA::RangeSegment(0, ncell), [=] RAJA_DEVICE (int cell) {
                 matrix_vector_code_kouter_atomics(cell, nlayers, lhs, x, ncell_3d, matrix,
                     ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
            });
        }
    }

    else if (memcmp(traverse,"colouring-kinner",12)==0){
        printf("Starting computation with colouring and kinner\n");
        for (int iter = 1; iter <= niters; iter ++){
            for (int colour = 0; colour < ncolour; colour ++){
                // for (int ccell = 0; ccell < ncp_colour[colour]; ccell ++){
                //     int cell = cmap[ccell*4+colour] - 1;
                //     matrix_vector_code_kinner(cell, nlayers, lhs, x, ncell_3d, matrix_kinner,
                //         ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
                // }
                RAJA::forall<policy>(RAJA::RangeSegment(0, ncp_colour[colour]), [=] RAJA_DEVICE (int ccell) {
                    int cell = cmap[ccell*4+colour] - 1;
                    matrix_vector_code_kinner(cell, nlayers, lhs, x, ncell_3d, matrix_kinner,
                        ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
                });
            }
        }
    }
    else if (memcmp(traverse,"colouring",6)==0){
        printf("Colouring traversing version\n");
        for (int iter = 1; iter <= niters; iter ++){
            for (int colour=0; colour < ncolour; colour ++){
                // for (int ccell = 0; ccell < ncp_colour[colour]; ccell ++){
                //     int cell = cmap[ccell*4+colour] - 1;
                //     matrix_vector_code_kouter(cell, nlayers, lhs, x, ncell_3d, matrix,
                //         ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
                // }
                RAJA::forall<policy>(RAJA::RangeSegment(0, ncp_colour[colour]), [=] RAJA_DEVICE (int ccell) {
                    int cell = cmap[ccell*4+colour] - 1;
                    matrix_vector_code_kouter(cell, nlayers, lhs, x, ncell_3d, matrix,
                        ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
                });
            }
        }
    }
}
