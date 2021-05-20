# Environment variables for OpenCL for Xilinx FPGAs

export FPGA_PLATFORM=xilinx_u200_xdma_201830_2
export FPGA_EXECUTION_TARGET=hw # Options: sw_emu | hw_emu | hw
export OCL_DEVICE_FLAGS="-O2 --report_level 1 --config xilinx.cfg"
export OCL_COMPILER="v++ --compile"
export OCL_LINKER="v++ --link"

# Host compiler flags
export OPENCL_LIBS="-lOpenCL"
export OPENCL_INCLUDE=""


# Configurable parameters
echo "Set FPGA_EXECUTION_TARGET to switch the target platform:"
echo "    sw_emu - for software emulation"
echo "    hw_emu - for hardware emulation"
echo "    hw - for FPGA execution. (default)"
