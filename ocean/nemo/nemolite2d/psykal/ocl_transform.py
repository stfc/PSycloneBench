def trans(psy):
    ''' Python script intended to be passed to PSyclone's generate()
    function via the -s option. Applies OpenCL to the Schedule so
    that PSyclone will generate an OpenCL PSy layer. '''

    from psyclone.psyGen import TransInfo
    tinfo = TransInfo()
    globaltrans = tinfo.get_trans_name('KernelGlobalsToArguments')
    cltrans = tinfo.get_trans_name('OCLTrans')

    # Get the invoke
    schedule = psy.invokes.get('invoke_0').schedule

    # Remove the globals from inside each kernel
    for idx, kern in enumerate(schedule.kernels()):
        globaltrans.apply(kern)

    # Transform invoke to OpenCL
    cltrans.apply(schedule)

    return psy
