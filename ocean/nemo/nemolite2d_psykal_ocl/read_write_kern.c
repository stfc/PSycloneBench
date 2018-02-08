
#ifdef __OPENCL_VERSION__

#ifdef INTELFPGA_CL
#pragma OPENCL EXTENSION cl_intel_channels : enable
channel float ssh_channel __attribute__((depth(10)));
#endif

__kernel void channel_write(int nx, int ny,
			    __global double* restrict ssha){
    for(int i=0; i<10; i++){
      write_channel_intel(ssh_channel, (float)(ssha[i]));
    }
    //mem_fence(CLK_CHANNEL_MEM_FENCE);
}


__kernel void channel_read(int nx, int ny,
			   __global double* restrict sshn){
  for(int i=0; i<1; i++){
    sshn[i] = (double)read_channel_intel(ssh_channel);
  }
  //mem_fence(CLK_CHANNEL_MEM_FENCE);
}

#endif
