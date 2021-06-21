''' Python script intended to be passed to PSyclone's generate()
function via the -s option. Applies OpenCL to the Schedule so
that PSyclone will generate an OpenCL PSy layer. '''

import os
from psyclone.psyGen import TransInfo
from psyclone.domain.gocean.transformations import \
    GOMoveIterationBoundariesInsideKernelTrans
from psyclone.configuration import Config


# Global variables to configure the PSyclone OpenCL generation:

# If FUCTIONAL_PARALLELISM is enabled the kernels will be dispatched to
# multiple command queues and run asynchronously when possible, otherwise
# all kernels will be run in-order in command queue number 1.
FUCTIONAL_PARALLELISM = True

# If MOVE_BOUNDARIES is True (mandatory for now) the start and stop boundaries
# of each loop will be moved from the PSy-layer to inside the kernel.
MOVE_BOUNDARIES = True

# If XILINX_CONFIG_FILE is True, this script will also generate a Xiling .cfg
# configuration file in the output directory.
XILINX_CONFIG_FILE = False

# The TILING parameter sets the number of kernel iterations that will be run
# together by a single kernel execution.
TILING = 64

def trans(psy):
    ''' Transform the schedule for OpenCL generation '''

    # Import transformations
    tinfo = TransInfo()
    globaltrans = tinfo.get_trans_name('KernelGlobalsToArguments')
    move_boundaries_trans = GOMoveIterationBoundariesInsideKernelTrans()
    cltrans = tinfo.get_trans_name('OCLTrans')

    # Get the invoke routine
    schedule = psy.invokes.get('invoke_0').schedule

    # Map the kernels by their name to different OpenCL queues. The multiple
    # command queues can be executed concurrently while each command queue
    # executes in-order its kernels. This provides functional parallelism
    # when kernels don't have dependencies between them.
    qmap = {
        'continuity_code': 1,
        'momentum_u_code': 2,
        'momentum_v_code': 3,
        'bc_ssh_code': 1,
        'bc_solid_u_code': 2,
        'bc_solid_v_code': 3,
        'bc_flather_u_code': 2,
        'bc_flather_v_code': 3,
        'field_copy_code': 1,
        'next_sshu_code': 1,
        'next_sshv_code': 1
        }

    # Remove global variables from inside each kernel, pass the boundary
    # values as arguments to the kernel and set the OpenCL work size to 64,
    # which is required for performance (with OpenCL < 1.2 this requires
    # the resulting application to be executed with DL_ESM_ALIGNMENT=64)
    for kern in schedule.kernels():
        print(kern.name)
        globaltrans.apply(kern)
        if MOVE_BOUNDARIES:
            move_boundaries_trans.apply(kern)
        if FUCTIONAL_PARALLELISM:
            kern.set_opencl_options({'local_size': TILING,
                                     'queue_number': qmap[kern.name]})
        else:
            kern.set_opencl_options({'local_size': TILING})

    # Transform invoke to OpenCL
    cltrans.apply(schedule)

    if XILINX_CONFIG_FILE:
        # Create a Xilinx Compiler Configuration file
        path = Config.get().kernel_output_dir
        with open(os.path.join(path, "xilinx.cfg"), "w") as cfgfile:
            cfgfile.write("# Xilinx FPGA configuration file\n")
            # cfgfile.write("[connectivity]\n")
            # cfgfile.write("# Create 2 CU of the given kernels\n")
            # cfgfile.write("nk=continuity_code:2\n")
            # cfgfile.write("nk=momentum_u_code:2\n")
            # cfgfile.write("nk=momentum_v_code:2\n")
            # cfgfile.write("nk=bc_ssh_code:2\n")

            # cfgfile.write("\n[hls]\n")
            # cfgfile.write("# Assign CUs to different SLRs\n")
            # cfgfile.write("slr=momentum_u_code_1:SLR0\n")
            # cfgfile.write("slr=momentum_u_code_2:SLR0\n")
            # cfgfile.write("slr=momentum_v_code_1:SLR2\n")
            # cfgfile.write("slr=momentum_v_code_2:SLR2\n")

    return psy
