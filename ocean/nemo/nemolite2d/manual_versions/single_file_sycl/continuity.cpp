#include "kernels.hpp"

template<typename TYPE>
void Continuity<TYPE>::operator()(item<1> item){
    auto i = item.get_id();
    ht[i] = hu[i] + hv[i]; 
}
