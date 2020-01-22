def trans(psy):
    ''' Python script intended to be passed to PSyclone's generate()
    function via the -s option. Applies OpenCL to the Schedule so
    that PSyclone will generate an OpenCL PSy layer. '''

    from psyclone.psyGen import TransInfo
    tinfo = TransInfo()
    globaltrans = tinfo.get_trans_name('KernelGlobalsToArguments')
    cltrans = tinfo.get_trans_name('OCLTrans')
    verbose = False

    # First invoke
    schedule = psy.invokes.get('invoke_0').schedule

    # Map of OpenCL options
    ocl_options =  {'local_size': 4}

    # Remove the globals from inside each kernel
    for idx, kern in enumerate(schedule.kernels()):
        print(kern.name)
        globaltrans.apply(kern)

    cltrans.apply(schedule)


    return psy
    # Attach options object to Kernel Call and Kernel Schedule
    for idx, kern in enumerate(schedule.kernels()):
        kern.set_opencl_options(ocl_options)
        kschedule = kern.get_kernel_schedule()
        if verbose:
            print("")
            kschedule.symbol_table.view()
            kschedule.view()
        kschedule.opencl_options = ocl_options

    return psy
