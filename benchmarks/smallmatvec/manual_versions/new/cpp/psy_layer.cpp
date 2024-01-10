#include <cstdio>
#include <cstdlib>
#include <cstring>


void matrix_vector_code_optimised(int cell, int nlayers, double *lhs, double *x, int ncell_3d, double *matrix, int ndf1, int undf1,
		int *map1, int ndf2, int undf2, int *map2){

	int ik = (cell-1)*nlayers;
	for (int df = 1; df <= ndf1; df ++){
		for (int df2 =1; df2 <= ndf2; df2 ++){
			for (int k = 1; k <= nlayers; k ++){
			//	printf("%d %d %d %d \n", df, df2, k, ik);
			//	printf("%d \n", (ik+k-1)*ndf1*ndf2+((df2-1)*ndf1)+df-1);
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
		for (int iter=1; iter <= niters; iter ++){
			for (int cell=1; cell <= ncell; cell ++){
				matrix_vector_code_optimised(cell, nlayers, lhs, x, ncell_3d, matrix,
						ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
			}
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
        else if (memcmp(traverse,"colouring",6)==0){
		printf("Colouring traversing version\n");
		for (int iter = 1; iter <= niters; iter ++){
			for (int colour=1; colour <= ncolour; colour ++){
				for (int ccell = 1; ccell <= ncp_colour[colour-1]; ccell ++){
					int cell = cmap[(ccell-1)*4+(colour-1)];
					matrix_vector_code_optimised(cell, nlayers, lhs, x, ncell_3d, matrix,
							ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);

				}
			}
		}
	}

}
