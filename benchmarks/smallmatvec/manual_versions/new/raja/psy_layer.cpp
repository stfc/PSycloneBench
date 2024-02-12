#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "RAJA/RAJA.hpp"

RAJA_DEVICE void matrix_vector_code_optimised(int cell, int nlayers, double *lhs, double *x, int ncell_3d, double *matrix, int ndf1, int undf1,
		int *map1, int ndf2, int undf2, int *map2){

	int ik = (cell-1)*nlayers;
	for (int df = 1; df <= ndf1; df ++){
		for (int df2 =1; df2 <= ndf2; df2 ++){
			for (int k = 1; k <= nlayers; k ++){
				lhs[map1[df-1]+k-2] = lhs[map1[df-1]+k-2] + matrix[(ik+k-1)*ndf1*ndf2+((df2-1)*ndf1)+df-1] * x[map2[df2-1]+k-2];
			}

		}
	}
}


void matrix_vector_code_kinner_atomics(int cell, int nlayers, double *lhs, double *x, int ncell_3d,
		double *matrix, int ndf1, int undf1, int *map1, int ndf2, int undf2, int *map2){

	double *x_e = (double*)malloc(sizeof(double)*ndf2);
	double *lhs_e = (double*)malloc(sizeof(double)*ndf1);
	int m2;
	int m1;

	for (int df2 = 1; df2 <= ndf2; df2 ++){
		m2 = map2[df2-1];
		for (int df = 1; df <= ndf1; df ++){
			m1 = map1[df-1];
			// #pragma omp simd
			for (int k=1; k <= nlayers; k++){
				lhs[m1+k-2]= lhs[m1+k-2] + matrix[(k-1) + (cell-1)*nlayers*ndf1*ndf2 + (df2-1)*nlayers*ndf1 + (df-1)*nlayers] * x[m2+k-2];
			}
		}

	}
}

extern "C" void c_psy_layer(char *traverse, int niters, int ncell, int nlayers,
		int ncell_3d, double *lhs, int *map_lhs,
		int ndf_lhs, int undf_lhs, double *matrix, double *matrix_kinner,
		double *x, int *map_x, int ndf_x, int undf_x,
		int ncolour, int *ncp_colour, int *cmap) {


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

// Define RAJA typed variables for the loop index variables 
// RAJA_INDEX_VALUE_T iter
// RAJA_INDEX_VALUE_T cell

	if (memcmp(traverse,"linear-kinner",11)==0){
		printf("Starting computation with linear and kinner\n");
		for (int iter = 1; iter <= niters; iter ++){
			for (int cell=1; cell <= ncell; cell ++){
				matrix_vector_code_kinner_atomics(cell, nlayers, lhs, x, ncell_3d, matrix_kinner,
						ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
			}


		}
	}
	else if (memcmp(traverse,"linear",5)==0){
		printf("Linear traversing Version\n");
		//using EXEC_POL =
		 // RAJA::KernelPolicy<
		   //   RAJA::statement::For<0,policy,
		     //   RAJA::statement::Lambda<0>
		       // >
	             // >;
	        for (int iter=1; iter <= niters; iter ++){
		//	RAJA::kernel<EXEC_POL>(RAJA::make_tuple(RAJA::TypedRangeSegment<int>(1, ncell)),
                  //      [=] (int cell){
			 //for (int cell=1; cell <= ncell; cell ++){
		        RAJA::forall<policy>(RAJA::RangeSegment(1, ncell), [=] RAJA_DEVICE (int cell) {
				matrix_vector_code_optimised(cell, nlayers, lhs, x, ncell_3d, matrix,
						ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
			}//);
		//}//);
			);
		}
	}

	else if (memcmp(traverse,"colouring-kinner",12)==0){
                printf("Starting computation with colouring and kinner\n");
                for (int iter = 1; iter <= niters; iter ++){
                         for (int colour=1; colour <= ncolour; colour ++){
                                 for (int ccell = 1; ccell <= ncp_colour[colour-1]; ccell ++){
                                          int cell = cmap[(ccell-1)*4+(colour-1)];
                                          matrix_vector_code_kinner_atomics(cell, nlayers, lhs, x, ncell_3d, matrix_kinner,
                                                ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
                                 }
                        }
                }
        }

// RAJA_INDEX_VALUE_T(KIDX, int, "KIDX");
// RAJA_INDEX_VALUE_T(JIDX, int, "JIDX"); 
// RAJA_INDEX_VALUE_T(IIDX, int, "IIDX");

// Define named loop index integer types used in the nested loops
// RAJA_INDEX_VALUE_T(KIDX, int, "iter");
// RAJA_INDEX_VALUE_T(JIDX, int, "colour");
// RAJA_INDEX_VALUE_T(IIDX, int, "ccell");

// Define min and max intervals for each loop index
// constexpr int imin = 1;
// constexpr int imax = ncp_colour[colour-1];
// constexpr int jmin = 1;
// constexpr int jmax = ncolour;
// constexpr int kmin = 1;
// constexpr int kmax = niters;

// Define corresponding typed range segments
// RAJA::TypedRangeSegment<KIDX> KRange(kmin, kmax);
// RAJA::TypedRangeSegment<JIDX> JRange(jmin, jmax);
// RAJA::TypedRangeSegment<IIDX> IRange(imin, imax);

// RAJA Loop - this is just loop reordering
// using KJI_EXECPOL = RAJA::KernelPolicy<
//                       RAJA::statement::For<2, RAJA::seq_exec,    // k
//                         RAJA::statement::For<1, RAJA::seq_exec,  // j
//                           RAJA::statement::For<0, RAJA::seq_exec,// i
//                             RAJA::statement::Lambda<0>
//                           >
//                         >
//                       >
//                     >;
//  RAJA::kernel<KJI_EXECPOL>( RAJA::make_tuple(IRange, JRange, KRange),
//  [=] (IIDX i, JIDX j, KIDX k) {
//     printf( " (%d, %d, %d) \n", (int)(*i), (int)(*j), (int)(*k));
//  });

        else if (memcmp(traverse,"colouring",6)==0){
		printf("Colouring traversing version\n");
		// for (int iter = 1; iter <= niters; iter ++){
		// 	for (int colour=1; colour <= ncolour; colour ++){
		// 		for (int ccell = 1; ccell <= ncp_colour[colour-1]; ccell ++){
		// 			int cell = cmap[(ccell-1)*4+(colour-1)];
		// 			matrix_vector_code_optimised(cell, nlayers, lhs, x, ncell_3d, matrix,
		// 					ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);

		// 		}
		// 	}
		// }
	}

}
