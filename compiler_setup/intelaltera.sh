# Environment variables for OpenCL for Intel Altera FPGA compilation

export FPGA_BOARD=nalla_pcie

# Host compiler flags
export OPENCL_LIBS="-L${INTELFPGAOCLSDKROOT}/board/${FPGA_BOARD}/linux64/lib \
    -L${INTELFPGAOCLSDKROOT}/host/linux64/lib -Wl,--no-as-needed -lalteracl \
    -l${FPGA_BOARD}_mmd -lelf"
export OPENCL_INCLUDE="-I${INTELFPGAOCLSDKROOT}/host/include"

# Configurable parameters
export CL_CONTEXT_EMULATOR_DEVICE_ALTERA=1
echo "Switch on/off the CL_CONTEXT_EMULATOR_DEVICE_ALTERA envvar to enable or \
    disable the FPGA emulation"
