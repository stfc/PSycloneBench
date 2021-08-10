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

    # Setup the dependency tracking
    next_dependencies = []
    for child in schedule.children:
        next_dependencies.append(0)

    # Fill the next_dependencies array with the forward dependencies
    i = 0
    for child in schedule.children:
        next_dependencies[i] = child.forward_dependence()
        i = i + 1

    # Loop through the dependencies, convert them from nodes to indices
    i = 0
    for child in schedule.children:
        j = 0
        for child in schedule.children:
            if child == next_dependencies[i]:
                next_dependencies[i] = j
                break
            j = j + 1
        i = i + 1

    # Loop through dependencies, and remove unneccessary ones.
    # This is first done by looping through each dependence.
    # For each dependence, we then loop over all other dependences
    # between the nodes involved in that dependence. If this dependence
    # fulfills a following dependence, the following dependence is removed.
    # If this dependence is fulfilled by a following dependence, this dependence
    # is removed.
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

    # Now we have computed (what we think is) a minimal set of 
    # dependencies, add the taskwait directives required. We do
    # this in reverse order to ensure we don't invalidate the
    # positions in the schedule.
    for i in range(len(next_dependencies)-1, -1, -1):
        if next_dependencies[i] is not None:
            schedule.addchild(OMPTaskwaitDirective(), next_dependencies[i])


    single_trans = OMPSingleTrans()
    parallel_trans = OMPParallelTrans()
    single_trans.apply(schedule.children)
    parallel_trans.apply(schedule.children)


