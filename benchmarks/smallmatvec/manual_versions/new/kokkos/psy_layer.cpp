#include <cstdio>

#include <Kokkos_Core.hpp>

extern "C" void c_psy_layer(char *traverse, int niters, int ncell, int nlayers,
                            int ncell_3d, double *lhs, int *map_lhs,
                            int ndf_lhs, int undf_lhs, double *matrix,
                            double *x, int *map_x, int ndf_x, int undf_x,
                            int ncolour, int *ncp_colour, int *cmap) {

  printf("Kokkos Version\n");
  Kokkos::initialize();
  Kokkos::finalize();
}
