def trans(psy):
    ''' Python script intended to be passed to PSyclone's generate()
    function via the -s option. Applies OpenCL to the Schedule so
    that PSyclone will generate an OpenCL PSy layer. '''

    from psyclone.psyGen import TransInfo
    tinfo = TransInfo()
    cltrans = tinfo.get_trans_name('OCLTrans')

    # First invoke
    schedule = psy.invokes.get('invoke_0').schedule
    cltrans.apply(schedule)

    # Map of OpenCL options
    ocl_options =  {'local_size': 8}

    # Attach options object to Kernel Call and Kernel Schedule
    for idx, kern in enumerate(schedule.kernels()):
        kern.set_opencl_options(ocl_options)
        kschedule = kern.get_kernel_schedule()
        kschedule.opencl_options = ocl_options

    return psy
