#include <CL/sycl.hpp>

#include <array>
#include <iostream>
#include <CL/sycl/intel/fpga_extensions.hpp>

using namespace cl::sycl;

constexpr access::mode sycl_read = access::mode::read;
constexpr access::mode sycl_write = access::mode::write;
constexpr access::mode sycl_read_write = access::mode::read_write;

using real = cl_double;

class Kernel1;

static const size_t ARRAY_SIZE = 127 * 127;

/* Function Headers */
void initialize_field(std::array<real, ARRAY_SIZE>& arr, const real value);
void initialize_tmask(std::array<cl_int, ARRAY_SIZE>& tmask, int nx, int ny,
                     int xstart, int xstop, int ystart, int ystop);

void add_vectors_parallel(std::array<cl_int, ARRAY_SIZE>& sum_array,
        const std::array<cl_int, ARRAY_SIZE>& addend_array_1,
        const std::array<cl_int, ARRAY_SIZE>& addend_array_2);

/* NemoLite2D application */
int main() {

    /* NemoLite2D Parameters */
    int nx = 127;
    int ny = 127;
    int xstart = 1;
    int xstop = nx - 1;
    int ystart = 1;
    int ystop = ny - 1;
    int istep;
    int nsteps = 100;
    int ji, jj, ikern;
    int buff_size;
    real dep_const = 100.0;
    real dx = 1000.0, dy=1000.0;
    real rdt = 20.0;
    real visc = 0.1;
    real cbfr = 0.00015;

    /* Sea-surface height */
    std::array<real, ARRAY_SIZE> ssha, ssha_u, ssha_v, sshn, sshn_u, sshn_v;
    std::array<real, ARRAY_SIZE> hu, hv, ht, un, vn, ua, va;
    std::array<real, ARRAY_SIZE> gphiu, gphiv;
    std::array<real, ARRAY_SIZE> e1u, e1v, e1t, e2u, e2v, e2t, e12u, e12v, e12t;

    /* T-point mask */
    std::array<cl_int, ARRAY_SIZE> tmask;

    // Initialize array with a given value
    initialize_field(hu, dep_const);
    initialize_field(hv, dep_const);
    initialize_field(ht, dep_const);
    initialize_field(un, 0.01);
    initialize_field(vn, 0.0);
    initialize_field(sshn_u, cos((360.0*ji)/(real)nx));
    initialize_field(sshn_v, cos((360.0*ji)/(real)nx));
    initialize_field(sshn, sin((360.0*ji)/(real)nx));
    initialize_field(ssha, 0.0);
    initialize_field(e1u, dx);
    initialize_field(e1v, dx);
    initialize_field(e1t, dx);
    initialize_field(e2u, dy);
    initialize_field(e2v, dy);
    initialize_field(e2t, dy);
    initialize_field(e12u, dx*dy);
    initialize_field(e12v, dx*dy);
    initialize_field(e12t, dx*dy);
    initialize_field(gphiu, 50.0);
    initialize_field(gphiv, 50.0);

    // Select SYCL device selector
    #if defined(FPGA_EMU)
        intel::fpga_emulator_selector device_selector;
    #elif defined(FPGA)
        intel::fpga_selector device_selector;
	#elif defined(GPU)
        gpu_selector device_selector;
    #elif defined(CPU)
        cpu_selector device_selector;
    #else
        throw -1;
	#endif
    std::unique_ptr<queue> device_queue;

    // Catch device selector runtime error 
    try {
        device_queue.reset( new queue(device_selector) );
    } catch (cl::sycl::exception const& e) {
        std::cout << "Caught a synchronous SYCL exception:" << std::endl << e.what() << std::endl;
        std::cout << "If you are targeting an FPGA hardware, please "
                    "ensure that your system is plugged to an FPGA board that is set up correctly"
                    "and compile with -DFPGA" << std::endl;
        std::cout << "If you are targeting the FPGA emulator, compile with -DFPGA_EMULATOR." << std::endl;
        return -1;
    }
    
    // Print selected device information    
    std::cout << "Device: "
            << device_queue->get_device().get_info<info::device::name>()
            << std::endl;

    // The size of amount of memory that will be given to the buffer 
    range<1> num_items{ARRAY_SIZE};

    // Create SYCL buffers: data shared between the host and the device
    buffer<real, 1> hu_buf(hu.data(), num_items);

    // Submit work to the device queue
    device_queue->submit([&](handler& cgh) {
        // Define the accessors
        auto hu_accessor = hu_buf.template get_access<sycl_read_write>(cgh);

        // Use parallel_for to execute each work_item in parallel
        cgh.parallel_for<class Kernel1>(num_items, [=](id<1> wiID) {
            hu_accessor[wiID] = hu_accessor[wiID] + hu_accessor[wiID];
            });
    });

    std::cout << "Finished successfully" << std::endl;
    return 0;
}

void initialize_field(std::array<real, ARRAY_SIZE>& arr, const real value){
    for (size_t i = 0; i < ARRAY_SIZE; i++){
        arr[i] = value;
    }
}

void initialize_tmask(
        std::array<cl_int, ARRAY_SIZE>& tmask, int nx, int ny,
        int xstart, int xstop, int ystart, int ystop){
    
    int jj, ji, idx;

    for(jj=0;jj<ny;jj++){
        for(ji=0;ji<nx;ji++){
            tmask[jj*nx + ji] = 1;
        }
    }

    for(jj=0;jj<ny;jj++){
        idx = jj*nx;
        // West solid boundary
        for(ji=0; ji<xstart; ji++){
            tmask[idx+ji] = 0;
        }
        // East solid boundary
        for(ji=xstop; ji<nx; ji++){
            tmask[idx+ji] = 0;
        }
    }
    // Southern open boundary
    for(jj=0; jj<ystart; jj++){
        idx = jj*nx;
        for(ji=0;ji<nx;ji++){
            tmask[idx + ji] = -1;
        }
    }
    // North solid boundary
    for(jj=ystop; jj<ny; jj++){
        idx = jj*nx;
        for(ji=0;ji<nx;ji++){
            tmask[idx + ji] = 0;
        }
    }
}
