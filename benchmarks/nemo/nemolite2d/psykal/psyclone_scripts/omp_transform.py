''' Python script intended to be passed to PSyclone's generate()
function via the -s option. It applies OpenMP to every loop and
inlines all kernels in the schedule.'''

from psyclone.configuration import Config
from psyclone.domain.common.transformations import KernelModuleInlineTrans
from psyclone.psyGen import TransInfo
from psyclone.psyir.nodes import Loop, Routine


def trans(psy):
    ''' Transformation entry point '''
    config = Config.get()
    tinfo = TransInfo()
    parallel_loop_trans = tinfo.get_trans_name('GOceanOMPParallelLoopTrans')
    loop_trans = tinfo.get_trans_name('GOceanOMPLoopTrans')
    parallel_trans = tinfo.get_trans_name('OMPParallelTrans')
    module_inline_trans = KernelModuleInlineTrans()

    schedule = psy.walk(Routine)[0]

    # Inline all kernels in this Schedule
    for kernel in schedule.kernels():
        module_inline_trans.apply(kernel)

    # Apply the OpenMPLoop transformation to every child in the schedule or
    # OpenMPParallelLoop to every Loop if it has distributed memory.
    for child in schedule.children:
        # We need to ignore dependencies on '{u,v}a' because PSyclone correctly
        # spots that there is a dependence in the bc_flather_{u,v} kernel.
        # However, we know that practically this isn't a problem
        # because these boundary-condition kernels only update values
        # outside the domain.
        options = {}
        if child.kernels()[0].name == "bc_flather_v_code":
            options["ignore_dependencies_for"] = ["va%data"]
        if child.kernels()[0].name == "bc_flather_u_code":
            options["ignore_dependencies_for"] = ["ua%data"]
        if config.distributed_memory:
            if isinstance(child, Loop):
                parallel_loop_trans.apply(child, options=options)
        else:
            loop_trans.apply(child, options=options)

    if not config.distributed_memory:
        # If it is not distributed memory, enclose all of these loops
        # within a single OpenMP PARALLEL region
        parallel_trans.apply(schedule.children)

    return psy
