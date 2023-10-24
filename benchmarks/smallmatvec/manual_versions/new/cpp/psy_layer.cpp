#include <cstdio>
#include <cstdlib>
#include <cstring>

void matmul(double* lhs, double* matrix, int dim1, int dim2, double* vector){
        for (int i = 0; i <= dim1 - 1; i ++){
		lhs[i] = 0;
		for (int j = 0; j <= dim2 - 1; j ++){
		lhs[i] += matrix[j+i*dim2] * vector[j];
		}
	}	
}

void matrix_vector_code_original_atomic(int cell, int nlayers, double *lhs, double *x, int ncell_3d,
		double *matrix, int ndf1, int undf1, int *map1, int ndf2, int undf2, int *map2){

	double *x_e = (double*)malloc(sizeof(double)*ndf2);
	double *lhs_e = (double*)malloc(sizeof(double)*ndf1);
	int ik;



	for (int k = 0; k <= nlayers-1; k ++){
		for (int df =1; df <= ndf2; df ++){

			x_e[df-1] = x[map2[df-1]+k-1];			
		}
		ik = (cell-1)*nlayers + k + 1;
	        matmul(lhs_e, &matrix[ik* ndf1 * ndf2 - 1], ndf1, ndf2, x_e);	
		for (int df = 1; df <= ndf1; df ++){
			lhs[map1[df-1]+k-1] = lhs[map1[df-1]+k-1] + lhs_e[df-1];
		}

	}
}

extern "C" void c_psy_layer(char *traverse, int niters, int ncell, int nlayers,
		int ncell_3d, double *lhs, int *map_lhs,
		int ndf_lhs, int undf_lhs, double *matrix,
		double *x, int *map_x, int ndf_x, int undf_x,
		int ncolour, int *ncp_colour, int *cmap) {

	printf("CPP Version\n");
	if (memcmp(traverse,"linear",5)==0){
		printf("Lineal traversing Version\n");
		for (int iter=1; iter <= niters; iter ++){
			printf("%d\n", iter);
			for (int cell=1; cell <= ncell; cell ++){
				matrix_vector_code_original_atomic( cell, nlayers, lhs, x, ncell_3d, matrix,
						ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
			}
		}
	}



}
