''' Python script intended to be passed to PSyclone's generate() function
via the -s option. This script module-inline all kernels in the PSy-layer.'''

from psyclone.psyGen import TransInfo
from psyclone.psyir.transformations import LoopTiling2DTrans
from psyclone.psyir.nodes import Loop
from psyclone.domain.gocean.transformations import \
    GOMoveIterationBoundariesInsideKernelTrans

# If MOVE_BOUNDARIES is True the start and stop boundaries of each loop will
# be moved from the PSy-layer to a mask inside the kernel. This is useful in
# combination with tiling as then all loop trip counts are equal an all loops
# can be made exactly divisible by a single tiling size.
MOVE_BOUNDARIES = True

# The TILING parameter sets the number of kernel iterations that will be run
# together by a single kernel execution.
TILING = 1


def trans(psy):
    ''' Transformation script entry function '''

    tinfo = TransInfo()
    itrans = tinfo.get_trans_name('KernelModuleInline')

    schedule = psy.invokes.get('invoke_0').schedule
    move_boundaries_trans = GOMoveIterationBoundariesInsideKernelTrans()

    # Module-Inline all coded kernels in this Schedule
    for kernel in schedule.coded_kernels():
        itrans.apply(kernel)

        if MOVE_BOUNDARIES:
            move_boundaries_trans.apply(kernel)

    # Add tiling to all loops
    if TILING > 1:
        for loop in schedule.walk(Loop):
            if loop.loop_type == "outer":
                LoopTiling2DTrans().apply(
                    loop, {"tilesize": TILING,
                           "strategy": 'exaclty_divisible'})

    return psy
