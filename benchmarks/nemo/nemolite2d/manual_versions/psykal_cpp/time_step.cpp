#include <iostream>
#include <chrono>
#include <vector>
#include <cstdlib>


extern "C" void c_invoke_time_step(
        double * c_ssha_t,
        double * c_sshn_t,
        double * c_sshn_u,
        double * c_sshn_v,
        double * c_hu,
        double * c_hv,
        double * c_hn,
        double * c_vn,
        double * c_ua,
        double * c_ht,
        double * c_ssha_u,
        double * c_va,
        double * c_ssha_v,
        int istp,
        int nx,
        int ny){

#ifdef DEBUG
    std::cout << "Hello from C++ in iteration: " << istp << std::endl;
    std::cout << "Size of the arrays: " << nx << " " << ny << std::endl;
    auto start = std::chrono::system_clock::now();
#endif

    int nsize = nx; // Simplify to squares and ignore ny for now

    // Currently using std::vector with a raw c array to vector transformation
    // that requires an expensive allocation and copy operations.
    // A better no-copy solution is needed!
    std::vector<std::vector<double> > ssha_t(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > sshn_t(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > sshn_u(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > sshn_v(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > hu(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > hv(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > hn(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > vn(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > ua(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > ht(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > ssha_u(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > va(nsize, std::vector<double>(nsize));
    std::vector<std::vector<double> > ssha_v(nsize, std::vector<double>(nsize));

    for(int i=0; i < nsize; i++){
        for(int j=0; j < nsize; j++){
            ssha_t[i][j] = c_ssha_t[i * nsize + j];
            sshn_t[i][j] = c_sshn_t[i * nsize + j];
            sshn_u[i][j] = c_sshn_u[i * nsize + j];
            sshn_v[i][j] = c_sshn_v[i * nsize + j];
            hu[i][j] = c_hu[i * nsize + j];
            hv[i][j] = c_hv[i * nsize + j];
            hn[i][j] = c_hn[i * nsize + j];
            vn[i][j] = c_vn[i * nsize + j];
            ua[i][j] = c_ua[i * nsize + j];
            ht[i][j] = c_ht[i * nsize + j];
            ssha_u[i][j] = c_ssha_u[i * nsize + j];
            va[i][j] = c_va[i * nsize + j];
            ssha_v[i][j] = c_ssha_v[i * nsize + j];
        }
    }

        

#ifdef DEBUG
    auto stop = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<
                    std::chrono::milliseconds>(stop - start);
    std::cout << "Time to allocate and copy raw pointers to vectors: ";
    std::cout << elapsed.count() << " ms"<< std::endl;
    start = std::chrono::system_clock::now();
#endif

    for(int i=0; i < nsize; i++){
        for(int j=0; j < nsize; j++){
            ssha_t[i][j] += 1.0;
            sshn_t[i][j] += 2.0;
            sshn_u[i][j] += 3.0;
            sshn_v[i][j] += 4.0;
            hu[i][j] += 5.0;
            hv[i][j] += 6.0;
            hn[i][j] += 7.0;
            vn[i][j] += 8.0;
            ua[i][j] += 9.0;
            ht[i][j] += 10.0;
            ssha_u[i][j] += 11.0;
            va[i][j] += 12.0;
            ssha_v[i][j] += 13.0;
        }
    }

#ifdef DEBUG
    stop = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<
                    std::chrono::milliseconds>(stop - start);
    std::cout << "Time to run the time step: ";
    std::cout << elapsed.count() << " ms"<< std::endl;
#endif

    for(int i=0; i < nsize; i++){
        for(int j=0; j < nsize; j++){
            c_ssha_t[i * nsize + j] = ssha_t[i][j];
            c_sshn_t[i * nsize + j] = sshn_t[i][j];
            c_sshn_u[i * nsize + j] = sshn_u[i][j];
            c_sshn_v[i * nsize + j] = sshn_v[i][j];
            c_hu[i * nsize + j] = hu[i][j];
            c_hv[i * nsize + j] = hv[i][j];
            c_hn[i * nsize + j] = hn[i][j];
            c_vn[i * nsize + j] = vn[i][j];
            c_ua[i * nsize + j] = ua[i][j];
            c_ht[i * nsize + j] = ht[i][j];
            c_ssha_u[i * nsize + j] = ssha_u[i][j];
            c_va[i * nsize + j] = va[i][j];
            c_ssha_v[i * nsize + j] = ssha_v[i][j];
        }
    }

#ifdef DEBUG
    stop = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<
                    std::chrono::milliseconds>(stop - start);
    std::cout << "Time to copy back the vectors to raw pointers: ";
    std::cout << elapsed.count() << " ms"<< std::endl;
#endif


}

/*
int main (int argc, char *argv[]) {

    int nsize = 1000;
    std::cout << "Quick test of c_invoke_time_step with size: ";
    std::cout << nsize << std::endl;

    // Allocate C style arrays to emulate Fortran call
    double * ssha_t =  (double*) malloc(nsize*nsize*sizeof(double));
    double * sshn_t =  (double*) malloc(nsize*nsize*sizeof(double));
    double * sshn_u =  (double*) malloc(nsize*nsize*sizeof(double));
    double * sshn_v =  (double*) malloc(nsize*nsize*sizeof(double));
    double * hu =  (double*) malloc(nsize*nsize*sizeof(double));
    double * hv =  (double*) malloc(nsize*nsize*sizeof(double));
    double * hn =  (double*) malloc(nsize*nsize*sizeof(double));
    double * vn =  (double*) malloc(nsize*nsize*sizeof(double));
    double * ua =  (double*) malloc(nsize*nsize*sizeof(double));
    double * ht =  (double*) malloc(nsize*nsize*sizeof(double));
    double * ssha_u =  (double*) malloc(nsize*nsize*sizeof(double));
    double * va =  (double*) malloc(nsize*nsize*sizeof(double));
    double * ssha_v =  (double*) malloc(nsize*nsize*sizeof(double));
    int istp = 0;

    c_invoke_time_step(ssha_t, sshn_t, sshn_u, sshn_v, hu, hv, hn,
                       vn, ua, ht, ssha_u, va, ssha_v, istp, nsize);

}*/
