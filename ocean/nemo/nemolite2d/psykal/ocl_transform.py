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
    ocl_options = [
        {'local_size':'8', 'queue_number':'2'},
        {'local_size':'8', 'queue_number':'3'},
        {'local_size':'8', 'queue_number':'4'},
        {'local_size':'8', 'queue_number':'2'},
        {'local_size':'8', 'queue_number':'3'},
        {'local_size':'8', 'queue_number':'4'},
        {'local_size':'8', 'queue_number':'3'},
        {'local_size':'8', 'queue_number':'4'},
        ]

    # Attach options object to Kernel Call and Kernel Schedule
    for idx, kern in enumerate(schedule.kernels()):
        print("Kernel " + str(idx) + " - " + kern.name + " : "
              + str(ocl_options[idx]))
        kern.opencl_options = ocl_options[idx]
        kschedule = kern.get_kernel_schedule()
        kschedule.opencl_options = ocl_options[idx]

    # Second invoke
    schedule = psy.invokes.get('invoke_1').schedule
    cltrans.apply(schedule)

    ocl_options = [
        {'local_size':'8', 'queue_number':'2'},
        {'local_size':'8', 'queue_number':'3'},
        {'local_size':'8', 'queue_number':'4'},
        {'local_size':'8', 'queue_number':'4'},
        {'local_size':'8', 'queue_number':'4'},
        ]

    for idx, kern in enumerate(schedule.kernels()):
        print("Kernel " + str(idx) + " - " + kern.name + " : "
              + str(ocl_options[idx]))
        kern.opencl_options = ocl_options[idx]
        kschedule = kern.get_kernel_schedule()
        kschedule.opencl_options = ocl_options[idx]

    return psy
