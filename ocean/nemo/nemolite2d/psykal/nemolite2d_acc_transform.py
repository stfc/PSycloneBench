
'''Python script intended to be passed to PSyclone's generate()
funcation via the -s option. Performs OpenACC transformations. '''


def trans(psy):
    ''' Take the supplied psy object, apply OpenACC transformations
    to the schedule of invoke_0 and return the new psy object '''
    from psyclone.transformations import OpenACCParallelTrans, \
        OpenACCDataTrans
    atrans = OpenACCParallelTrans()
    dtrans = OpenACCDataTrans()

    invoke = psy.invokes.get('invoke_0')
    schedule = invoke.schedule
    # schedule.view()

    # Apply the OpenMP Loop transformation to *every* loop
    # in the schedule
    from psyclone.psyGen import Loop
    for child in schedule.children:
        if isinstance(child, Loop):
            newschedule, _ = atrans.apply(child)
            schedule = newschedule

    newschedule, _ = dtrans.apply(schedule)

    invoke.schedule = newschedule
    newschedule.view()
    return psy


if __name__ == "__main__":
    from psyclone.parse import parse
    from psyclone.psyGen import PSyFactory
    API = "gocean1.0"
    FILENAME = "nemolite2d_alg.f90"
    _, INVOKEINFO = parse(FILENAME,
                          api=API,
                          invoke_name="invoke")
    PSY = PSyFactory(API).create(INVOKEINFO)
    # print PSY.invokes.names

    NEW_PSY = trans(PSY)

    print NEW_PSY.gen
