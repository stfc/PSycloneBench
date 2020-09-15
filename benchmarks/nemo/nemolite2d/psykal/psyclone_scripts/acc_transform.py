
'''Python script intended to be passed to PSyclone's generate()
funcation via the -s option. Performs OpenACC transformations. '''


def trans(psy):
    ''' Take the supplied psy object, apply OpenACC transformations
    to the schedule of invoke_0 and return the new psy object '''
    from psyclone.psyGen import TransInfo
    tinfo = TransInfo()
    parallel_trans = tinfo.get_trans_name('ACCParallelTrans')
    loop_trans = tinfo.get_trans_name('ACCLoopTrans')
    enter_data_trans = tinfo.get_trans_name('ACCEnterDataTrans')
    routine_trans = tinfo.get_trans_name('ACCRoutineTrans')
    glo2arg_trans = tinfo.get_trans_name('KernelGlobalsToArguments')
    inline_trans = tinfo.get_trans_name('KernelModuleInline')

    invoke = psy.invokes.get('invoke_0')
    schedule = invoke.schedule

    # Apply the OpenACC Loop transformation to *every* loop
    # in the schedule
    from psyclone.psyir.nodes import Loop
    for child in schedule.children:
        if isinstance(child, Loop):
            newschedule, _ = loop_trans.apply(child, {"collapse": 2})
            schedule = newschedule

    # Put all of the loops in a single parallel region
    parallel_trans.apply(schedule)

    # Add an enter-data directive
    enter_data_trans.apply(schedule)

    # Apply ACCRoutineTrans to each kernel, which it requires to remove
    # any global variables first.
    for kern in schedule.coded_kernels():
        glo2arg_trans.apply(kern)
        routine_trans.apply(kern)
        inline_trans.apply(kern)

    # schedule.view()
    return psy
