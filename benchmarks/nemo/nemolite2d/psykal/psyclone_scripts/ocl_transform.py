''' Python script intended to be passed to PSyclone's generate()
function via the -s option. Applies OpenCL to the Schedule so
that PSyclone will generate an OpenCL PSy layer. '''

from psyclone.psyGen import TransInfo

def trans(psy):
    ''' Transform the schedule for OpenCL generation '''
    tinfo = TransInfo()
    globaltrans = tinfo.get_trans_name('KernelGlobalsToArguments')
    cltrans = tinfo.get_trans_name('OCLTrans')

    # Get the invoke
    schedule = psy.invokes.get('invoke_0').schedule

    # Remove global variables from inside each kernel and set the OpenCL
    # work size to 64, which is required for performance (currently
    # this can cause memory issues with arbitrary problem sizes).
    for kern in schedule.kernels():
        globaltrans.apply(kern)
        kern.set_opencl_options({'local_size': 64})

    # Transform invoke to OpenCL
    cltrans.apply(schedule)

    return psy
