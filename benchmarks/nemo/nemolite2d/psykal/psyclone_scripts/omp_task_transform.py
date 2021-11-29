''' Python script intended to be passed to PSyclone's generate()
function via the -s option. It applies OpenMP to every loop and
inlines all kernels in the schedule.'''

from psyclone.psyGen import TransInfo
from psyclone.psyir.nodes import Loop
from psyclone.configuration import Config
from psyclone.transformations import OMPParallelTrans, OMPSingleTrans, \
                                     OMPTaskloopTrans, KernelModuleInlineTrans
from psyclone.psyir.transformations import OMPTaskwaitTrans
from psyclone.psyir.nodes import OMPTaskloopDirective, OMPTaskwaitDirective, \
                                 OMPDirective, OMPParallelDirective
from psyclone.psyGen import HaloExchange


def trans(psy):
    '''Transformation entry point'''
    config = Config.get()

    schedule = psy.invokes.get('invoke_0').schedule



    loop_trans = OMPTaskloopTrans(grainsize=32, nogroup=True)
    wait_trans = OMPTaskwaitTrans()

    module_inline_trans = KernelModuleInlineTrans()


    # Inline all kernels in this Schedule
    for kernel in schedule.kernels():
        module_inline_trans.apply(kernel)

    for child in schedule.children:
        if isinstance(child, Loop):
            loop_trans.apply(child)

    single_trans = OMPSingleTrans()
    parallel_trans = OMPParallelTrans()
    sets = []
    if config.distributed_memory:
        # Find all of the groupings of taskloop and taskwait directives. Each of these
        # groups needs its own parallel+single regions.
        next_start = 0
        next_end = 0
        for index in range(len(schedule.children)):
            child = schedule.children[index]
            if isinstance(child, OMPTaskloopDirective) or isinstance(child, OMPTaskwaitDirective):
                next_end = next_end + 1
            elif not isinstance(child, OMPDirective):
                if next_start == index:
                    next_end = index + 1
                    next_start = index + 1
                else:
                    sets.append((next_start, next_end))
                    next_end = index + 1
                    next_start = index+1
            else:
                next_end = next_end + 1
        if next_start <= index:
            sets.append((next_start, index+1))
        sets.reverse()
        for next_set in sets:
            single_trans.apply(schedule[next_set[0]:next_set[1]])
            parallel_trans.apply(schedule[next_set[0]])
        for child in schedule.children:
            if isinstance(child, OMPParallelDirective):
                wait_trans.apply(child)
    else:
        single_trans.apply(schedule.children)
        parallel_trans.apply(schedule.children)
        wait_trans.apply(schedule.children[0])
