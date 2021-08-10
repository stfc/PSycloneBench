''' Python script intended to be passed to PSyclone's generate()
function via the -s option. It applies OpenMP to every loop and
inlines all kernels in the schedule.'''

from psyclone.psyGen import TransInfo
from psyclone.psyir.nodes import Loop
from psyclone.configuration import Config
from psyclone.transformations import OMPParallelTrans, OMPSingleTrans
from psyclone.transformations import OMPTaskloopTrans, KernelModuleInlineTrans
from psyclone.psyir.nodes import OMPTaskwaitDirective


def trans(psy):
    '''Transformation entry point'''
    config = Config.get()

    schedule = psy.invokes.get('invoke_0').schedule

    if config.distributed_memory:
        print("Distributed memory not yet implemented")
        quit()

    loop_trans = OMPTaskloopTrans(grainsize=32)

    module_inline_trans = KernelModuleInlineTrans()

    # Inline all kernels in this Schedule
    for kernel in schedule.kernels():
        module_inline_trans.apply(kernel)

    next_dependencies = []
    for child in schedule.children:
        next_dependencies.append(0)

    i = 0

    for child in schedule.children:
        next_dependencies[i] = child.forward_dependence()
        i = i + 1

    i = 0
    for child in schedule.children:
        j = 0
        for child in schedule.children:
            if child == next_dependencies[i]:
                next_dependencies[i] = j
                break
            j = j + 1
        i = i + 1

    for i in range(len(next_dependencies)):
        if next_dependencies[i] is not None:
            next_dependence = next_dependencies[i]
            for j in range(i+1, next_dependence):
                if j >= next_dependence:
                    break
                if next_dependencies[j] is not None:
                    if next_dependencies[j] < next_dependence:
                        next_dependence = next_dependencies[j]
                        next_dependencies[i] = None
                    if next_dependencies[j] > next_dependence:
                        next_dependencies[j] = None

    for child in schedule.children:
        loop_trans.apply(child)
    print(next_dependencies)
    # TODO: Sort out dependencies!
    for i in range(len(next_dependencies)-1, -1, -1):
        if next_dependencies[i] is not None:
            schedule.addchild(OMPTaskwaitDirective(), next_dependencies[i])


    single_trans = OMPSingleTrans()
    parallel_trans = OMPParallelTrans()
    single_trans.apply(schedule.children)
    parallel_trans.apply(schedule.children)


