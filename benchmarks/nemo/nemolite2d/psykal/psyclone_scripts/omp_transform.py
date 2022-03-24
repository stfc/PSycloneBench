''' Python script intended to be passed to PSyclone's generate()
function via the -s option. It applies OpenMP to every loop and
inlines all kernels in the schedule.'''

from psyclone.psyGen import TransInfo
from psyclone.psyir.nodes import Loop
from psyclone.configuration import Config
from psyclone.domain.gocean.transformations import \
    GOMoveIterationBoundariesInsideKernelTrans
from psyclone.psyir.transformations import LoopTiling2DTrans

# If MOVE_BOUNDARIES is True the start and stop boundaries of each loop will
# be moved from the PSy-layer to inside the kernel.
MOVE_BOUNDARIES = True

# The TILING parameter sets the number of kernel iterations that will be run
# together by a single kernel execution.
TILING = 1


def trans(psy):
    ''' Transformation entry point '''
    config = Config.get()
    tinfo = TransInfo()
    parallel_loop_trans = tinfo.get_trans_name('GOceanOMPParallelLoopTrans')
    loop_trans = tinfo.get_trans_name('GOceanOMPLoopTrans')
    parallel_trans = tinfo.get_trans_name('OMPParallelTrans')
    module_inline_trans = tinfo.get_trans_name('KernelModuleInline')
    move_boundaries_trans = GOMoveIterationBoundariesInsideKernelTrans()

    schedule = psy.invokes.get('invoke_0').schedule

    # Inline all kernels in this Schedule
    for kernel in schedule.kernels():
        module_inline_trans.apply(kernel)

        if MOVE_BOUNDARIES:
            move_boundaries_trans.apply(kernel)

    if TILING > 1:
        for loop in schedule.walk(Loop):
            if loop.loop_type == "outer":
                LoopTiling2DTrans().apply(loop,
                                          {"tilesize": TILING,
                                           "strategy": 'exactly_divisible'})

    # Apply the OpenMPLoop transformation to every child in the schedule or
    # OpenMPParallelLoop to every Loop if it has distributed memory.
    for loop in schedule.walk(Loop, stop_type=(Loop)):
        if loop.loop_type == "outer":

            if config.distributed_memory:
                parallel_loop_trans.apply(loop)
            else:
                loop_trans.apply(loop)

    if not config.distributed_memory:
        # If it is not distributed memory, enclose all of these loops
        # within a single OpenMP PARALLEL region
        parallel_trans.apply(schedule.children)

    return psy
