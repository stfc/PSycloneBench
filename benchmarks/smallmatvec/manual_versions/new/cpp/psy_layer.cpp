#include <cstdio>
#include <cstdlib>
#include <cstring>


void matrix_vector_code_kouter(
        int cell, int nlayers, double *restrict lhs, double *restrict x, int ncell_3d, double *restrict matrix, int ndf1, int undf1,
        int *restrict map1, int ndf2, int undf2, int *restrict map2
){

    int ik = cell*nlayers;
    for (int df2 = 0; df2 < ndf2; df2 ++){
        int m2 = map2[df2] - 1; // -1 because map2 contains fortran 1-indexing references
        for (int df = 0; df < ndf1; df ++){
            int m1 = map1[df] - 1; // -1 because map2 contains fortran 1-indexing references
            #pragma omp simd
            for (int k = 0; k < nlayers; k ++){
                lhs[m1+k] = lhs[m1+k] + matrix[(ik+k)*ndf1*ndf2 + df2*ndf1 + df] * x[m2+k];
            }
        }
    }
}

void matrix_vector_code_kouter_atomic(
        int cell, int nlayers, double *restrict lhs, double *restrict x, int ncell_3d, double *restrict matrix, int ndf1, int undf1,
        int *restrict map1, int ndf2, int undf2, int *restrict map2
){

    int ik = cell*nlayers;
    for (int df2 = 0; df2 < ndf2; df2 ++){
        int m2 = map2[df2] - 1; // -1 because map2 contains fortran 1-indexing references
        for (int df = 0; df < ndf1; df ++){
            int m1 = map1[df] - 1; // -1 because map2 contains fortran 1-indexing references
            for (int k = 0; k < nlayers; k ++){
                #pragma omp atomic
                lhs[m1+k] = lhs[m1+k] + matrix[(ik+k)*ndf1*ndf2 + df2*ndf1 + df] * x[m2+k];
            }
        }
    }
}

void matrix_vector_code_kinner(
        int cell, int nlayers, double *restrict lhs, double *restrict x, int ncell_3d,
        double *restrict matrix, int ndf1, int undf1, int *restrict map1, int ndf2, int undf2, int *restrict map2
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

void matrix_vector_code_kinner_atomics(
        int cell, int nlayers, double *restrict lhs, double *restrict x, int ncell_3d,
        double *restrict matrix, int ndf1, int undf1, int *restrict map1, int ndf2, int undf2, int *restrict map2
){

    for (int df2 = 0; df2 < ndf2; df2 ++){
        int m2 = map2[df2] - 1; // -1 because map2 contains fortran 1-indexing references
        for (int df = 0; df < ndf1; df ++){
            int m1 = map1[df] - 1; // -1 because map2 contains fortran 1-indexing references
            for (int k = 0; k < nlayers; k++){
                #pragma omp atomic
                lhs[m1+k]= lhs[m1+k] + matrix[k + cell*nlayers*ndf1*ndf2 + df2*nlayers*ndf1 + df*nlayers] * x[m2+k];
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
#ifdef TARGET_GPU                        
            #pragma omp target loop
#else
            #pragma omp parallel for
#endif
            for (int cell=0; cell < ncell; cell ++){
                matrix_vector_code_kinner_atomics(cell, nlayers, lhs, x, ncell_3d, matrix_kinner,
                        ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
            }


        }
    }
    else if (memcmp(traverse,"linear",5)==0){
        printf("Linear traversing Version\n");
        for (int iter=1; iter <= niters; iter ++){
#ifdef TARGET_GPU
            #pragma omp target loop
#else
            #pragma omp parallel for
#endif
            for (int cell=0; cell < ncell; cell ++){
                matrix_vector_code_kouter_atomic(cell, nlayers, lhs, x, ncell_3d, matrix,
                        ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
            }
        }
    }
    else if (memcmp(traverse,"colouring-kinner",12)==0){
        printf("Starting computation with colouring and kinner\n");
        for (int iter = 1; iter <= niters; iter ++){
            for (int colour = 0; colour < ncolour; colour ++){
#ifdef TARGET_GPU
                #pragma omp target loop
#else
                #pragma omp parallel for
#endif
                for (int ccell = 0; ccell < ncp_colour[colour]; ccell ++){
                    int cell = cmap[ccell*4 + colour] - 1;
                    matrix_vector_code_kinner(cell, nlayers, lhs, x, ncell_3d, matrix_kinner,
                            ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);
                }
            }
        }
    }
    else if (memcmp(traverse,"colouring",6)==0){
        printf("Colouring traversing version\n");
        for (int iter = 1; iter <= niters; iter ++){
            for (int colour = 0; colour < ncolour; colour ++){
#ifdef TARGET_GPU
                #pragma omp target loop
#else
                #pragma omp parallel for
#endif
                for (int ccell = 0; ccell < ncp_colour[colour]; ccell ++){
                    int cell = cmap[ccell*4 + colour] - 1;
                    matrix_vector_code_kouter(cell, nlayers, lhs, x, ncell_3d, matrix,
                            ndf_lhs, undf_lhs, &map_lhs[cell*ndf_lhs], ndf_x, undf_x, &map_x[cell*ndf_x]);

                }
            }
        }
    }
}
