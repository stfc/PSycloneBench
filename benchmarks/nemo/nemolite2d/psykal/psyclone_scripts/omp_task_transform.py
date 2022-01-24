''' Python script intended to be passed to PSyclone's generate()
function via the -s option. It applies OpenMP tasking to every loop
and inlines all kernels in the schedule.'''

from psyclone.psyir.nodes import Loop
from psyclone.configuration import Config
from psyclone.transformations import OMPParallelTrans, OMPSingleTrans, \
                                     OMPTaskloopTrans, KernelModuleInlineTrans
from psyclone.psyir.transformations import OMPTaskwaitTrans
from psyclone.psyir.nodes import OMPTaskloopDirective, OMPTaskwaitDirective, \
                                 OMPDirective, OMPParallelDirective


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
    if not config.distributed_memory:
        single_trans.apply(schedule.children)
        parallel_trans.apply(schedule.children)
        wait_trans.apply(schedule.children[0])
        return
    # Find all of the groupings of taskloop and taskwait directives. Each of
    # these groups needs its own parallel+single regions. This makes sure we
    # don't apply OpenMP transformations to the Halo Exchange operations.
    next_start = 0
    next_end = 0
    idx = 0
    for idx, child in enumerate(schedule.children):
        # Loop through through the schedule until we find a non-OpenMP
        # node, and extend the grouping until we do.
        if isinstance(child, (OMPTaskloopDirective, OMPTaskwaitDirective)):
            next_end = next_end + 1
        elif not isinstance(child, OMPDirective):
            # If we find a non OpenMP directive, if we're currently in a
            # grouping of OpenMP directives then we stop, and add it to
            # the set of groupings. Otherwise we just skip over this
            # node.
            if next_start == idx:
                next_end = idx + 1
                next_start = idx + 1
            else:
                sets.append((next_start, next_end))
                next_end = idx + 1
                next_start = idx + 1
        else:
            next_end = next_end + 1
    # If currently in a grouping of directives, add it to the list
    # of groupings
    if next_start <= idx:
        sets.append((next_start, idx+1))
    # Start from the last grouping to keep indexing correct,
    # so reverse the ordering
    sets.reverse()
    # For each of the groupings of OpenMP directives, surround them
    # with an OpenMP Single and an OpenMP Parallel directive set.
    for next_set in sets:
        single_trans.apply(schedule[next_set[0]:next_set[1]])
        parallel_trans.apply(schedule[next_set[0]])
    # Finally, we loop over the OMPParallelDirectives, and apply the
    # OMPTaskWaitTrans to ensure correctness.
    for child in schedule.children:
        if isinstance(child, OMPParallelDirective):
            wait_trans.apply(child)
