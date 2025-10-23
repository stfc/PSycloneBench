''' Python script intended to be passed to PSyclone's generate() function
via the -s option. This script module-inline all kernels in the PSy-layer.'''

from psyclone.domain.common.transformations import KernelModuleInlineTrans
from psyclone.psyir.nodes import Node, Routine


def trans(psy: Node):
    '''Entry point for PSyIR transformation. This script module-inlines
    every user-supplied kernel that is called.

    '''
    itrans = KernelModuleInlineTrans()

    schedule = psy.walk(Routine)[0]

    # Module-Inline all coded kernels in this Schedule
    for kernel in schedule.coded_kernels():
        itrans.apply(kernel)

    return psy
