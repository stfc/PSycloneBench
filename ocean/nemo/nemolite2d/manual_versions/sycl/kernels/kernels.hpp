#ifndef NEMOLITE_KERNELS_H
#define NEMOLITE_KERNELS_H

#include <CL/sycl.hpp>

using namespace cl::sycl;

constexpr access::mode sycl_read = access::mode::read;
constexpr access::mode sycl_write = access::mode::write;
constexpr access::mode sycl_read_write = access::mode::read_write;

template<typename TYPE>
using read_accessor = accessor<TYPE, 1, sycl_read, access::target::global_buffer>;
template<typename TYPE>
using write_accessor = accessor<TYPE, 1, sycl_write, access::target::global_buffer>;
template<typename TYPE>
using read_write_accessor = accessor<TYPE, 1, sycl_read_write, access::target::global_buffer>;



template<typename TYPE> class Continuity {
private:
    read_accessor<TYPE> hu;
    read_accessor<TYPE> hv;
    read_write_accessor<TYPE> ht;

public:
    Continuity(read_accessor<TYPE> set_hu,
               read_accessor<TYPE> set_hv,
               read_write_accessor<TYPE> set_ht):
      hu(set_hu), hv(set_hv), ht(set_ht) { };

    void operator()(item<1> item){
        auto i = item.get_id();
        ht[i] = 1;// hu[i] + hv[i]; 
    };
};

#endif
