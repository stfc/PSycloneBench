''' Python script intended to be passed to PSyclone's generate()
function via the -s option. '''

from psyclone.psyGen import TransInfo

def trans(psy):
    ''' Transformation script entry function '''

    tinfo = TransInfo()
    itrans = tinfo.get_trans_name('KernelModuleInline')

    schedule = psy.invokes.get('invoke_0').schedule

    # Inline all kernels in this Schedule
    for kernel in schedule.kernels():
        itrans.apply(kernel)

    return psy
