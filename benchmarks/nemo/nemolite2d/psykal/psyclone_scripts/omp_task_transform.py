''' Python script intended to be passed to PSyclone's generate()
function via the -s option. It applies OpenMP tasking to every loop
and inlines all kernels in the schedule.'''

from psyclone.psyir.nodes import Loop
from psyclone.configuration import Config
from psyclone.transformations import OMPParallelTrans, OMPSingleTrans
from psyclone.domain.common.transformations import KernelModuleInlineTrans
from psyclone.psyir.transformations.fold_conditional_return_expressions_trans \
        import FoldConditionalReturnExpressionsTrans
from psyclone.psyir.transformations import OMPTaskTrans, ChunkLoopTrans, InlineTrans
from psyclone.psyir.nodes import OMPParallelDirective, OMPTaskDirective, \
                                 OMPDirective, Call, IntrinsicCall
from psyclone.transformations import \
    KernelImportsToArguments


def trans(psy):
    '''Transformation entry point'''
    config = Config.get()

    schedule = psy.invokes.get('invoke_0').schedule
    schedule.symbol_table.resolve_imports()
    print(schedule.symbol_table.get_symbols()['ssha_t'].datatype)
    print(schedule.symbol_table.get_symbols()['ssha_t'].datatype.datatype)
    resolved = schedule.symbol_table.get_symbols()['ssha_t'].datatype.resolve_deferred()
    print(resolved.datatype)
    print("----------")

    loop_trans = ChunkLoopTrans()
    task_trans = OMPTaskTrans()
    inline_trans = InlineTrans()
#    imports_to_arguments = KernelImportsToArguments()
    cond_trans = FoldConditionalReturnExpressionsTrans()
    module_inline_trans = KernelModuleInlineTrans()

    # Inline all kernels in this Schedule
    for kernel in schedule.kernels():
#        imports_to_arguments.apply(kernel)
        module_inline_trans.apply(kernel)
        cond_trans.apply(kernel.get_kernel_schedule())
        kernel.lower_to_language_level()
    for call in schedule.walk(Call):
#        import pdb
#        pdb.set_trace()
        if isinstance(call, IntrinsicCall):
            continue
        inline_trans.apply(call)
        #pass

    applications = 0
    for child in schedule.children[:]:
        if isinstance(child, Loop):
            #print(child.debug_string())
            loop_trans.apply(child)
            
            assert isinstance(child.children[3].children[0], Loop)
            task_trans.apply(child, {"force": True})
            applications = applications + 1

    print(f"Applied task transformation {applications} times.")

    single_trans = OMPSingleTrans()
    parallel_trans = OMPParallelTrans()
    sets = []
    if not config.distributed_memory:
        single_trans.apply(schedule.children)
        parallel_trans.apply(schedule.children)
        print(f"Found {len(schedule.walk(OMPTaskDirective))} task directives.")
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
        if isinstance(child, (OMPTaskDirective)):
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
