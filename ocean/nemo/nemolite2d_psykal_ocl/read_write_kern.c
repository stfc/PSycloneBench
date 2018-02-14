
#pragma OPENCL EXTENSION cl_intel_channels : enable

channel int ssh_channel __attribute__((depth(10)));

__kernel void channel_write(int nx, int ny,
			    __global double* restrict ssha){
    for(int i=0; i<100000; i++){
      write_channel_intel(ssh_channel, i); //(float)(ssha[i]));
    }
}


__kernel void channel_read(int nx, int ny,
			   __global double* restrict sshn){
  int j;
  for(int i=0; i<100000; i++){
    j = read_channel_intel(ssh_channel);
  }
  sshn[0] = (double)j;
}
