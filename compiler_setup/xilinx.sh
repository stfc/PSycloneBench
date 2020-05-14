# Environment variables for OpenCL for Xilinx FPGAs

export PLATFORM = xilinx_u200_xdma_201830_2
export EXECUTION_TARGET = sw_emu # Options: sw_emu | hw_emu | hw
export OCL_COMPILER = "v++ --compile --platform ${PLATFORM} --target ${EXECUTION_TARGET}"
export OCL_LINKER = "v++ --link --platform ${PLATFORM} --target ${EXECUTION_TARGET}"

# Host compiler flags
export OPENCL_LIBS="-lOpenCL"
export OPENCL_INCLUDE=""


# Configurable parameters
echo "Use the EXECUTION_TARGET env var to switch the target platform:  \
    sw_emu - for software emulation, hw_emu - for hardware emulation, \
    hw - for FPGA execution."
