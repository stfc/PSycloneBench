#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "RAJA/RAJA.hpp"
#include "chai/ManagedArray.hpp"
#include "umpire/ResourceManager.hpp"
#include "umpire/strategy/QuickPool.hpp"

#if defined(USE_CUDA_POLICY)
    using policy = RAJA::cuda_exec<256, false>;
    using atomic = RAJA::cuda_atomic;
#elif defined(USE_OPENMP_POLICY)
    using policy = RAJA::omp_parallel_for_exec;   
    using atomic = RAJA::omp_atomic;
#else
    using policy = RAJA::loop_exec;
    using atomic = RAJA::omp_atomic;
#endif

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
		    RAJA::atomicAdd<atomic>(&lhs[m1+k], matrix[(ik+k)*ndf1*ndf2 + df2*ndf1 + df] * x[m2+k]);
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
		    RAJA::atomicAdd<atomic>(&lhs[m1+k], matrix[k + cell*nlayers*ndf1*ndf2 + df2*nlayers*ndf1 + df*nlayers] * x[m2+k]);
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
    printf("Using RAJA cuda version \n");
#elif defined(USE_OPENMP_POLICY)
    printf("Using RAJA OpenMP version \n");
#else
    printf("Using RAJA sequential version \n");
#endif

    auto& rm = umpire::ResourceManager::getInstance();
    auto pool = rm.makeAllocator<umpire::strategy::QuickPool>("pool", rm.getAllocator("UM"));
    double * lhs_um = static_cast<double*>(pool.allocate(undf_lhs * sizeof(double)));
    double * x_um = static_cast<double*>(pool.allocate(undf_x * sizeof(double)));
    int * map_lhs_um = static_cast<int*>(pool.allocate(ndf_lhs * ncell * sizeof(int)));
    int * map_x_um = static_cast<int*>(pool.allocate(ndf_x * ncell * sizeof(int)));
    double * matrix_um = static_cast<double*>(pool.allocate(ndf_lhs * ndf_x * ncell_3d * sizeof(double)));
    printf("Umpire Pool allocation complete \n");

    RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0, undf_lhs), [=] (int idx) {
        lhs_um[idx] = lhs[idx];
    });
    RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0, ndf_lhs * ncell), [=] (int idx) {
        map_lhs_um[idx] = map_lhs[idx];
    });
    printf("Copy from C array to UM lhs array complete \n");
    RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0, undf_x), [=] (int idx) {
        x_um[idx] = x[idx];
    });
    RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0, ndf_x * ncell), [=] (int idx) {
        map_x_um[idx] = map_x[idx];
    });
    printf("Copy from C array to UM x array complete \n");
    RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0, ndf_lhs * ndf_x * ncell_3d), [=] (int idx) {
        matrix_um[idx] = matrix[idx];
    });
    printf("Copy from C array to UM matrix array complete \n");

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
            // RAJA::forall<policy>(RAJA::RangeSegment(0, ncell), [=] RAJA_DEVICE (int cell) {
            //      matrix_vector_code_kouter_atomics(cell, nlayers, lhs_um, x_um, ncell_3d, matrix_um,
            //          ndf_lhs, undf_lhs, &map_lhs_um[cell*ndf_lhs], ndf_x, undf_x, &map_x_um[cell*ndf_x]);
            // });
            RAJA::launch<RAJA::LaunchPolicy<RAJA::cuda_launch_t<false>>>( // For CPU: omp_launch_t
                RAJA::LaunchParams(
                    RAJA::Teams(ncell),
                    RAJA::Threads(1,ndf_lhs,nlayers)),
                    [=] RAJA_HOST_DEVICE (RAJA::LaunchContext ctx) {
                        // for (int cell=0; cell < ncell; cell ++){
                        RAJA::loop<RAJA::LoopPolicy<RAJA::cuda_block_x_direct>>(ctx, RAJA::RangeSegment(0, ncell), // For CPU: omp_for
                            [&] (int cell){
                                int ik = cell*nlayers;
                                RAJA::loop<RAJA::LoopPolicy<RAJA::loop_exec>>(ctx, RAJA::RangeSegment(0, ndf_x), // For CPU: omp_for
                                [&] (int df2){
                                // for (int df2 = 0; df2 < ndf_x; df2 ++){ // cuda_thread_x_direct, for cpu: loop/sequential
                                    int m2 = map_x_um[cell*ndf_x + df2] - 1; // -1 because map2 contains fortran 1-indexing references
                                    RAJA::loop<RAJA::LoopPolicy<RAJA::cuda_thread_y_direct>>(ctx, RAJA::RangeSegment(0, ndf_lhs), // For CPU: omp_for
                                    [&] (int df){
                                    // for (int df = 0; df < ndf_lhs; df ++){ //cuda_thread_y_direct
                                        int m1 = map_lhs_um[cell*ndf_lhs + df] - 1; // -1 because map2 contains fortran 1-indexing references
                                        RAJA::loop<RAJA::LoopPolicy<RAJA::cuda_thread_z_direct>>(ctx, RAJA::RangeSegment(0, nlayers), // For CPU: omp_for
                                        [&] (int k){
                                        // for (int k = 0; k < nlayers; k ++){ // cuda_thread_z_direct
                                            RAJA::atomicAdd<atomic>(&lhs_um[m1+k], matrix_um[(ik+k)*ndf_lhs*ndf_x + df2*ndf_lhs + df] * x_um[m2+k]);
                                        });
                                    });
                                });
                            });
                    });
            // for (int cell=0; cell < ncell; cell ++){
            //     int ik = cell*nlayers;
            //     for (int df2 = 0; df2 < ndf_x; df2 ++){
            //         int m2 = map_x_um[df2] - 1; // -1 because map2 contains fortran 1-indexing references
            //         for (int df = 0; df < ndf_lhs; df ++){
            //             int m1 = map_lhs_um[df] - 1; // -1 because map2 contains fortran 1-indexing references
            //             for (int k = 0; k < nlayers; k ++){
            // 		    RAJA::atomicAdd<atomic>(&lhs_um[m1+k], matrix_um[(ik+k)*ndf_lhs*ndf_x + df2*ndf_lhs + df] * x_um[m2+k]);
            //             }
            //         }
            //     }
            // }

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
    RAJA::forall<RAJA::loop_exec>(RAJA::RangeSegment(0, undf_lhs), [=] (int idx) {
        lhs[idx] = lhs_um[idx];
    });
}
