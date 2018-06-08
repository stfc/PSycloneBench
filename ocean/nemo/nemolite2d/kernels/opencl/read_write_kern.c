
#pragma OPENCL EXTENSION cl_intel_channels : enable

channel int transfer __attribute__((depth(10)));

__kernel void channel_write(int nx, int ny,
			    __global int* restrict ssha){
  int j;
  for(int i=0; i<nx*ny; i++){
    //j = i;
    write_channel_intel(transfer, i); //(float)(ssha[i]));
  }
  //ssha[0] = j;
}


__kernel void channel_read(int nx, int ny,
			   __global int* restrict sshn){
  int j;
  for(int i=0; i<nx*ny; i++){
    //sshn[i] = i;
    sshn[i] = read_channel_intel(transfer);
  }
}
