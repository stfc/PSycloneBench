#include <iostream>
#include <vector>
#include <cstdlib>

// Kernels
#include "../../kernels/c_family/continuity_kern.c"
#include "../../kernels/c_family/momentum_u_kern.c"
#include "../../kernels/c_family/momentum_v_kern.c"
#include "../../kernels/c_family/boundary_conditions_kern.c"
#include "../../kernels/c_family/time_update_kern.c"

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
        int * tmask,
        double * area_t,
        double * area_u,
        double * area_v,
        double * dx_u,
        double * dx_v,
        double * dx_t,
        double * dy_u,
        double * dy_v,
        double * dy_t,
        double * gphiu,
        double * gphiv,
        // Scalars
        int istep,
        int internal_xstart,
        int internal_xstop,
        int internal_ystart,
        int internal_ystop,
        int width,
        double rdt,
        double cbfr,
        double visc,
        double omega,
        double d2r,
        double g
        ){

#define TASK_SIZE 32
//Each double loop over jj then ji is divided into a new structure.
//The outer loop, now loops over increments of 32 (initially chosen arbitrarily),
//and creates tasks of this size
//For the dependencies supplied to OmpSs, as opposed to putting the full memory for now,
//we compute the "task row"
//For the nth task for each jj+=32 loop, this is n.
//Tasks read to the task_rows appropriately (always either task_row, task_row-1, or task_row+1)
//Tasks only write to the task_row according to jj/TASK_SIZE
//TASK_SIZE is controlled by the define above
//We use the task_rows to control the dependencies between tasks.
//OmpSs/OpenMP task dependencies don't need to cover the actual memory addresses used.
//As long as all tasks use the same sized TASK_SIZE, the dependencies will be sound.
//We could write this to use the memory addresses directly, but this is more readable
//and could allow the dependence analysis to be easier, enabling the code to run faster

    // Continuity kernel (internal domain)
    // Try dividing continuity code into 32 rows at a time
    for(int jj = internal_ystart; jj <= internal_ystop; jj+=TASK_SIZE){
      int stop = std::min(jj+TASK_SIZE-1, internal_ystop);
      int task_row = jj/TASK_SIZE;
      #pragma omp task firstprivate(jj, task_row, stop) \
                output(ssha_t[task_row]) input(sshn_u[task_row], \
                                               hu[task_row], \
                                               un[task_row], \
                                               sshn_v[task_row], \
                                               sshn_v[task_row-1], \
                                               hv[task_row], \
                                               hv[task_row-1], \
                                               vn[task_row], \
                                               vn[task_row], \
                                               sshn_t[task_row])
                                                                    
      for( int task_j = jj; task_j <= stop; task_j++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            continuity_code(ji, task_j, width, ssha_t, sshn_t, sshn_u, sshn_v, \
                hu, hv, un, vn, rdt, area_t);
        }
      }
    }

    // Momentum_u kernel (internal domain u points)
    // Try dividing momentum_u_code into 32 rows at a time
    for(int jj = internal_ystart; jj <= internal_ystop; jj+=TASK_SIZE){
      int stop = std::min(jj+TASK_SIZE-1, internal_ystop);
      int task_row = jj/TASK_SIZE;
      #pragma omp task firstprivate(jj, task_row, stop) \
                       inout( ua[ task_row ] ) \
                          in( un[ task_row], \
                              un[ task_row-1], \
                              vn[ task_row ],  \
                              vn[ task_row + 1], \
                              vn[ task_row - 1], \
                              hu[ task_row ], \
                              hu[ task_row + 1], \
                              hu[ task_row - 1], \
                              hv[ task_row ], \
                              hv[ task_row + 1], \
                              hv[ task_row - 1], \
                              ht[ task_row ], \
                              ht[ task_row + 1], \
                              ssha_u[ task_row ], \
                              sshn_t[ task_row ], \
                              sshn_t[ task_row + 1 ], \
                              sshn_t[ task_row - 1 ], \
                              sshn_v[ task_row ], \
                              sshn_v[ task_row + 1 ], \
                              sshn_v[ task_row - 1 ] )
      for( int task_j =jj; task_j <= stop; task_j++){
        for(int ji = internal_xstart; ji <= internal_xstop - 1; ji++){
            momentum_u_code(ji, task_j, width, ua, un, vn, hu, hv, ht, ssha_u, \
                sshn_t, sshn_u, sshn_v, tmask, dx_u, dx_v, dx_t, dy_u, dy_t, \
                area_u, gphiu, rdt, cbfr, visc, omega, d2r, g);
        }
      }
    }

    // Momentum_v kernel (internal domain v points)
    for(int jj = internal_ystart; jj <= internal_ystop - 1; jj+=TASK_SIZE){
      int stop = std::min(jj+TASK_SIZE-1, internal_ystop-1);
      int task_row = jj/TASK_SIZE;
      #pragma omp task firstprivate(jj, task_row, stop) \
                  inout( va[task_row] ) \
                     in( un[task_row], \
                         un[task_row - 1], \
                         un[task_row + 1], \
                         vn[task_row], \
                         vn[task_row - 1], \
                         vn[task_row + 1], \
                         hu[task_row], \
                         hu[task_row - 1], \
                         hu[task_row + 1], \
                         hv[task_row], \
                         hv[task_row - 1], \
                         hv[task_row + 1], \
                         ht[task_row], \
                         ht[task_row + 1], \
                         ssha_v[task_row], \
                         sshn_t[task_row], \
                         sshn_t[task_row + 1], \
                         sshn_u[task_row], \
                         sshn_u[task_row-1], \
                         sshn_u[task_row+1], \
                         sshn_v[task_row], \
                         sshn_v[task_row-1], \
                         sshn_v[task_row+1])
      for( int task_j =jj; task_j <= stop; task_j++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            momentum_v_code(ji, task_j, width, va, un, vn, hu, hv, ht, ssha_v, \
                sshn_t, sshn_u, sshn_v, tmask, dx_v, dx_t, dy_u, dy_v, dy_t, \
                area_v, gphiv, rdt, cbfr, visc, omega, d2r, g);
        }
      }
    }

    // Boundary conditions bc_ssh kernel (internal domain)
    for(int jj = internal_ystart; jj <= internal_ystop; jj+=TASK_SIZE){
      int stop = std::min(jj+TASK_SIZE-1, internal_ystop);
      int task_row = jj/TASK_SIZE;
      #pragma omp task firstprivate(jj, task_row, stop) out( ssha_t[task_row] )
      for( int task_j =jj; task_j <= stop; task_j++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            bc_ssh_code(ji, task_j, width, istep, ssha_t, tmask, rdt);
        }
      }
    }

    // Boundary conditions bc_solid_u kernel (whole domain but top x boundary)
    for(int jj = internal_ystart - 1; jj <= internal_ystop + 1; jj+=TASK_SIZE){
      int stop = std::min(jj+TASK_SIZE-1, internal_ystop+1);
      int task_row = jj/TASK_SIZE;
      #pragma omp task firstprivate(jj, task_row, stop) out( ua[task_row] )
      for( int task_j =jj; task_j <= stop; task_j++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop; ji++){
            bc_solid_u_code(ji, task_j, width, ua, tmask);
        }
      }
    }

    // Boundary conditions bc_solid_v kernel (whole domain but top y boundary)
    for(int jj = internal_ystart - 1; jj <= internal_ystop; jj+=TASK_SIZE){
      int stop = std::min(jj+TASK_SIZE-1, internal_ystop);
      int task_row = jj/TASK_SIZE;
      #pragma omp task firstprivate(jj, task_row, stop) out( ua[task_row] )
      for( int task_j =jj; task_j <= stop; task_j++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop + 1; ji++){
            bc_solid_v_code(ji, task_j, width, va, tmask);
        }
      }
    }
    // Boundary conditions bc_flather_u kernel (whole domain but top x boundary)
    for(int jj = internal_ystart - 1; jj <= internal_ystop + 1; jj+=TASK_SIZE){
      int stop = std::min(jj+TASK_SIZE-1, internal_ystop);
      int task_row = jj/TASK_SIZE;
      #pragma omp task firstprivate(jj, task_row, stop) out( ua[task_row] ) \
                                                         in( ua[task_row] , \
                                                             ua[task_row-1], \
                                                             ua[task_row+1], \
                                                             hu[task_row], \
                                                             sshn_u[task_row], \
                                                             sshn_u[task_row+1], \
                                                             sshn_u[task_row-1] )
               
                                                                      
      for( int task_j =jj; task_j <= stop; task_j++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop; ji++){
            bc_flather_u_code(ji, task_j, width, ua, hu, sshn_u, tmask, g);
        }
      }
    }

//Since the following code is not in a task, we must wait for the previous work to complete
//so we need this taskwait.
#pragma omp taskwait
    // Boundary conditions bc_flather_v kernel (whole domain but top y boundary)
    // Can't be parallelised over j, so no task (messes up dependence analysis for following tasks
    for(int jj = internal_ystart - 1; jj <= internal_ystop; jj++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop + 1; ji++){
            bc_flather_v_code(ji, jj, width, va, hv, sshn_v, tmask, g);
        }
    }

    // Copy 'next' fields to 'current' fields (whole domain)
    for(int jj = internal_ystart - 1; jj < internal_ystop + 1; jj+=TASK_SIZE){
      int stop = std::min(jj+TASK_SIZE-1, internal_ystop+1);
      int task_row = jj/TASK_SIZE;
        #pragma omp task firstprivate(jj, task_row, stop) out(un[task_row], \
                             vn[task_row], \
                             sshn_t[task_row]) \
                          in(ua[task_row], \
                             va[task_row], \
                             sshn_t[task_row])
      for( int task_j =jj; task_j <= stop; task_j++){
        for(int ji = internal_xstart - 1; ji <= internal_xstop + 1; ji++){
            int idx = task_j * width + ji;
            un[idx] = ua[idx];
            vn[idx] = va[idx];
            sshn_t[idx] = ssha_t[idx];
        }
      }
    }

    // Time update kernel (internal domain u points)
    for(int jj = internal_ystart; jj <= internal_ystop; jj+=TASK_SIZE){
      int stop = std::min(jj+TASK_SIZE-1, internal_ystop);
      int task_row = jj/TASK_SIZE;
        #pragma omp task firstprivate(jj, task_row, stop) out(sshn_u[task_row]) \
                                                           in(sshn_t[task_row], \
                                                              sshn_t[task_row+1])                                                              
      for( int task_j =jj; task_j <= stop; task_j++){
        for(int ji = internal_xstart; ji <= internal_xstop - 1; ji++){
            next_sshu_code(ji, task_j, width, sshn_u, sshn_t, tmask, area_t, area_u);
        }
      }
    }

    // Time update kernel (internal domain v points)
    for(int jj = internal_ystart; jj <= internal_ystop - 1; jj+=TASK_SIZE){
      int stop = std::min(jj+TASK_SIZE-1, internal_ystop - 1);
      int task_row = jj/TASK_SIZE;
        #pragma omp task firstprivate(jj, task_row, stop) out(sshn_v[task_row]) \
                                                           in(sshn_t[task_row], \
                                                              sshn_t[task_row+1])                                                              
      for( int task_j =jj; task_j <= stop; task_j++){
        for(int ji = internal_xstart; ji <= internal_xstop; ji++){
            next_sshv_code(ji, task_j, width, sshn_v, sshn_t, tmask, area_t, area_v);
        }
      }
    }
}
