''' Python script intended to be passed to PSyclone's generate()
function via the -s option. Applies OpenCL to the Schedule so
that PSyclone will generate an OpenCL PSy layer. '''

from psyclone.psyGen import TransInfo
from psyclone.domain.gocean.transformations import \
    GOMoveIterationBoundariesInsideKernelTrans


def trans(psy):
    ''' Transform the schedule for OpenCL generation '''
    tinfo = TransInfo()
    globaltrans = tinfo.get_trans_name('KernelGlobalsToArguments')
    move_boundaries_trans = GOMoveIterationBoundariesInsideKernelTrans()
    cltrans = tinfo.get_trans_name('OCLTrans')

    # Get the invoke
    schedule = psy.invokes.get('invoke_0').schedule

    # Remove global variables from inside each kernel, pass the boundary
    # values as arguments to the kernel and set the OpenCL work size to 64,
    # which is required for performance (with OpenCL < 1.2 this requires
    # the resulting application to be executed with DL_ESM_ALIGNMENT=64)
    for kern in schedule.kernels():
        globaltrans.apply(kern)
        move_boundaries_trans.apply(kern)
        kern.set_opencl_options({'local_size': 64})

    # Transform invoke to OpenCL
    cltrans.apply(schedule)

    return psy
