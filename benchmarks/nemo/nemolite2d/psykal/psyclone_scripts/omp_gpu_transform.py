''' Python script intended to be passed to PSyclone's generate()
function via the -s option. It applies OpenMP to every loop and
inlines all kernels in the schedule.'''

from psyclone.domain.common.transformations import KernelModuleInlineTrans
from psyclone.psyir.nodes import Routine, Loop
from psyclone.psyir.transformations import OMPTargetTrans
from psyclone.psyir.transformations import OMPLoopTrans
from psyclone.psyir.transformations import TransformationError
from psyclone.psyir.transformations import \
    FoldConditionalReturnExpressionsTrans
from psyclone.transformations import \
    KernelImportsToArguments, OMPDeclareTargetTrans


def trans(psy):
    ''' Transformation entry point '''
    loop_trans = OMPLoopTrans()
    target_trans = OMPTargetTrans()
    module_inline_trans = KernelModuleInlineTrans()
    imports_to_arguments = KernelImportsToArguments()
    fold_trans = FoldConditionalReturnExpressionsTrans()
    omp_declare = OMPDeclareTargetTrans()

    schedule = psy.invokes.invoke_list[0].schedule

    # Inline all kernels in this Schedule
    for kernel in schedule.kernels():
        imports_to_arguments.apply(kernel)
        fold_trans.apply(kernel.get_kernel_schedule())
        module_inline_trans.apply(kernel)
        kernel.lower_to_language_level()

    loop_trans.omp_directive = "teamsdistributeparalleldo"
    loop_trans.omp_schedule = "none"

    for outer_loop in schedule.walk(Loop, stop_type=Loop):
        try:
            loop_trans.apply(outer_loop)
            outer_loop.parent.parent.collapse = 2
            target_trans.apply(outer_loop.parent.parent)
        except TransformationError:
            pass

    for kernel_subroutine in schedule.parent.walk(Routine):
        if not kernel_subroutine.name.startswith("invoke_"):
            omp_declare.apply(kernel_subroutine)

    return psy