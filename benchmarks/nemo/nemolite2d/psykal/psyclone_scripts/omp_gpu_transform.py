''' Python script intended to be passed to PSyclone's generate()
function via the -s option. It applies OpenMP to every loop and
inlines all kernels in the schedule.'''

from psyclone.domain.common.transformations import KernelModuleInlineTrans
from psyclone.psyir.transformations import OMPTargetTrans
from psyclone.psyir.transformations import OMPLoopTrans
from psyclone.psyir.transformations import TransformationError
from psyclone.psyir.nodes import Loop


def trans(psy):
    ''' Transformation entry point '''
    loop_trans = OMPLoopTrans()
    target_trans = OMPTargetTrans()
    module_inline_trans = KernelModuleInlineTrans()

    schedule = psy.invokes.get('invoke_0').schedule

    # Inline all kernels in this Schedule
    for kernel in schedule.kernels():
        module_inline_trans.apply(kernel)

    loop_trans.omp_directive = "teamsdistributeparalleldo"

    for outer_loop in schedule.walk(Loop, stop_type=Loop):
        try:
            loop_trans.apply(outer_loop)
            outer_loop.parent.parent.collapse = 2
            target_trans.apply(outer_loop.parent.parent)
        except TransformationError:
            pass

    return psy
