''' Python script intended to be passed to PSyclone's generate()
function via the -s option. It applies OpenMP tasking to every loop
and inlines all kernels in the schedule.'''

from psyclone.psyir.nodes import Loop
from psyclone.configuration import Config
from psyclone.transformations import KernelModuleInlineTrans
from psyclone.psyir.transformations import OtterTaskloopTrans,\
        OtterParallelTrans, OtterSynchroniseRegionTrans
from psyclone.psyir.nodes import OtterNode, OtterTaskNode, OtterParallelNode


def trans(psy):
    '''Transformation entry point'''
    config = Config.get()

    schedule = psy.invokes.get('invoke_0').schedule

    loop_trans = OtterTaskloopTrans()

    module_inline_trans = KernelModuleInlineTrans()

    # Inline all kernels in this Schedule
    for kernel in schedule.kernels():
        module_inline_trans.apply(kernel)

    for child in schedule.children:
        if isinstance(child, Loop):
            loop_trans.apply(child)

    parallel_trans = OtterParallelTrans()
    sync_trans = OtterSynchroniseRegionTrans()
    sets = []
    if not config.distributed_memory:
        parallel_trans.apply(schedule.children)
        sync_trans.apply(schedule.children[0])
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
        if isinstance(child, (OtterTaskNode)):
            next_end = next_end + 1
        elif not isinstance(child, OtterNode):
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
        parallel_trans.apply(schedule[next_set[0]])
    # Finall do synchronisation inside each parallel region
    for child in schedule.walk(OtterParallelNode):
        sync_trans.apply(child)
