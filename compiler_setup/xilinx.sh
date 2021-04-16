# Environment variables for OpenCL for Xilinx FPGAs

export PLATFORM=xilinx_u200_xdma_201830_2
export EXECUTION_TARGET=hw # Options: sw_emu | hw_emu | hw
export XILINX_FLAGS="--jobs 8 -O2 --report_level 1 --config xilinx.cfg"
export OCL_COMPILER="v++ --compile --platform ${PLATFORM} --target ${EXECUTION_TARGET} ${XILINX_FLAGS}"
export OCL_LINKER="v++ --link --platform ${PLATFORM} --target ${EXECUTION_TARGET} ${XILINX_FLAGS}"

# Host compiler flags
export OPENCL_LIBS="-lOpenCL"
export OPENCL_INCLUDE=""


# Configurable parameters
echo "Set EXECUTION_TARGET to switch the target platform:"
echo "    sw_emu - for software emulation"
echo "    hw_emu - for hardware emulation"
echo "    hw - for FPGA execution. (default)"
