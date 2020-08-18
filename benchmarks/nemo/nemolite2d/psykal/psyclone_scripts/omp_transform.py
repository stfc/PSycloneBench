def trans(psy):
    ''' Python script intended to be passed to PSyclone's generate()
    function via the -s option. Applies OpenMP to every loop before
    enclosing them all within a single OpenMP PARALLEL region. '''

    from psyclone.psyGen import TransInfo
    tinfo = TransInfo()
    ltrans = tinfo.get_trans_name('GOceanOMPLoopTrans')
    rtrans = tinfo.get_trans_name('OMPParallelTrans')
    itrans = tinfo.get_trans_name('KernelModuleInline')

    schedule = psy.invokes.get('invoke_0').schedule

    # Inline all kernels in this Schedule
    for kernel in schedule.kernels():
        itrans.apply(kernel)

    # Apply the OpenMP Loop transformation to *every* loop
    # in the schedule
    for child in schedule.children:
        newschedule, _ = ltrans.apply(child)
        schedule = newschedule

    # Enclose all of these loops within a single OpenMP
    # PARALLEL region
    newschedule, _ = rtrans.apply(schedule.children)

    return psy
