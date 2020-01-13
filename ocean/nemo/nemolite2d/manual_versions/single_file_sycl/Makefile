CXX = dpcpp  # source /opt/intel/inteloneapi/setvars.sh
CXXFLAGS =
LDFLAGS = -lOpenCL -lsycl
SOURCES = nemolite2d.cpp


.PHONY: nemolite_cpu nemolite_gpu nemolite_fpga nemolite_fpga_emu run_emu clean 
all: nemolite_cpu

nemolite_cpu:
	$(CXX) -DCPU $(CXXFLAGS) $(SOURCES) $(LDFLAGS) -o nemolite_cpu.exe

nemolite_gpu:
	$(CXX) -DGPU $(CXXFLAGS) $(SOURCES) $(LDFLAGS) -o nemolite_gpu.exe

nemolite_fpga:
	$(CXX) -DFPGA $(CXXFLAGS) -fintelfpga $(SOURCES) -o nemolite_fpga.exe -Xshardware

nemolite_fpga_emu:
	$(CXX) -DFPGA_EMU $(CXXFLAGS) -fintelfpga $(SOURCES) -o nemolite_fpga_emu.exe

clean: 
	rm -rf nemolite_cpu.exe nemolite_gpu.exe nemolite_fpga.exe nemolite_fpga_emu.exe
