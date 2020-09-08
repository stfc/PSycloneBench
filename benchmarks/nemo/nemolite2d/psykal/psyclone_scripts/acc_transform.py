
'''Python script intended to be passed to PSyclone's generate()
funcation via the -s option. Performs OpenACC transformations. '''


def trans(psy):
    ''' Take the supplied psy object, apply OpenACC transformations
    to the schedule of invoke_0 and return the new psy object '''
    from psyclone.psyGen import TransInfo
    tinfo = TransInfo()
    atrans = tinfo.get_trans_name('ACCParallelTrans')
    dtrans = tinfo.get_trans_name('ACCDataTrans')
    itrans = tinfo.get_trans_name('KernelModuleInline')

    invoke = psy.invokes.get('invoke_0')
    schedule = invoke.schedule
    # schedule.view()

    # Apply the OpenACC Loop transformation to *every* loop
    # in the schedule
    from psyclone.psyir.nodes import Loop
    for child in schedule.children:
        if isinstance(child, Loop):
            newschedule, _ = atrans.apply(child)
            schedule = newschedule

    dtrans.apply(schedule)

    schedule.view()
    return psy
