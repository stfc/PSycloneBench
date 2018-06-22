def trans(psy):
    ''' Python script intended to be passed to PSyclone's generate()
    function via the -s option. Applies OpenCL to the Schedule so
    that PSyclone will generate an OpenCL PSy layer. '''

    from psyclone.psyGen import TransInfo
    tinfo = TransInfo()
    cltrans = tinfo.get_trans_name('OCLTrans')

    schedule = psy.invokes.get('invoke_0').schedule
    # schedule.view()
    cltrans.apply(schedule)

    psy.invokes.get('invoke_0').schedule = newschedule
    return psy
