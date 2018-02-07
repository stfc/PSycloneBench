
#ifdef __OPENCL_VERSION__

#ifdef INTELFPGA_CL
#pragma OPENCL EXTENSION cl_intel_channels : enable
channel double ssh_channel __attribute__((depth(10)));
#endif

__kernel void channel_write(int nx, int ny,
			    __global double* restrict ssha){

    write_channel_intel(ssh_channel, ssha[0]);
    //mem_fence(CLK_CHANNEL_MEM_FENCE);
}


__kernel void channel_read(int nx, int ny,
			   __global double* restrict sshn){
  bool valid;
  sshn[0] = read_channel_intel(ssh_channel);//, &valid);
  //mem_fence(CLK_CHANNEL_MEM_FENCE);
}

#endif
