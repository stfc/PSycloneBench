
#ifdef __OPENCL_VERSION__

#ifdef INTELFPGA_CL
#pragma OPENCL EXTENSION cl_intel_channels : enable
channel double ssh_channel __attribute__((depth(100)));
#endif

__kernel void channel_write(int nx, int ny,
			    __global double* restrict ssha){

    write_channel_nb_intel(ssh_channel, ssha[0]);
}


__kernel void channel_read(int nx, int ny,
			   __global double* restrict sshn){
  bool valid;
  sshn[0] = read_channel_nb_intel(ssh_channel, &valid);
}

#endif
