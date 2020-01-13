
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// Physical Parameters
__constant static const double PI = 3.1415926535897932;
__constant static const double g = 9.80665 ;
__constant static const double omega = 7.292116e-05 ;
__constant static const double d2r = PI/180.0 ;

// Model parameters (NEEDS TO MATCH WITH NAMELIST)
__constant static const double rdt = 20.0 ;
__constant static const double cbfr = 0.00015 ;
__constant static const double visc = 0.1 ;

