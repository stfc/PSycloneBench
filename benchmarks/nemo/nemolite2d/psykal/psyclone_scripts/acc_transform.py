'''Python script intended to be passed to PSyclone's generate()
function via the -s option. Performs OpenACC transformations. '''

from psyclone.domain.common.transformations import KernelModuleInlineTrans
from psyclone.psyGen import TransInfo
from psyclone.psyir.nodes import Loop, Routine
from psyclone.transformations import (
    ACCEnterDataTrans, ACCLoopTrans, ACCParallelTrans, ACCRoutineTrans,
    KernelImportsToArguments)


def trans(psy):
    ''' Take the supplied psy object, apply OpenACC transformations
    to the schedule of invoke_0 and return the new psy object '''
    tinfo = TransInfo()
    parallel_trans = tinfo.get_trans_name('ACCParallelTrans')
    loop_trans = tinfo.get_trans_name('ACCLoopTrans')
    enter_data_trans = ACCEnterDataTrans()
    routine_trans = ACCRoutineTrans()
    glo2arg_trans = KernelImportsToArguments()
    inline_trans = KernelModuleInlineTrans()

    schedule = psy.walk(Routine)[0]

    # Apply the OpenACC Loop transformation to *every* loop
    # in the schedule
    for child in schedule.children:
        if isinstance(child, Loop):
            # We need to ignore dependencies on 'va' because PSyclone correctly
            # spots that there is a dependence in one of the boundary-condition
            # kernels. However, we know that practically this isn't a problem
            # because of the way the domain (mask) is configured.
            loop_trans.apply(child, {"collapse": 2,
                                     "ignore_dependencies_for": ["va"]})

    # Put all of the loops in a single parallel region
    parallel_trans.apply(schedule)

    # Add an enter-data directive
    enter_data_trans.apply(schedule)

    # Apply ACCRoutineTrans to each kernel, which also requires that
    # any global variables must be removed first.
    for kern in schedule.coded_kernels():
        glo2arg_trans.apply(kern)
        routine_trans.apply(kern)
        inline_trans.apply(kern)

    return psy
