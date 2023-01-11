''' Python script intended to be passed to PSyclone's generate() function
via the -s option. This script module-inline all kernels in the PSy-layer.'''

from psyclone.domain.common.transformations import KernelModuleInlineTrans


def trans(psy):
    ''' Transformation script entry function.

    :param psy: The PSy layer object to apply transformations to.
    :type psy: :py:class:`psyclone.psyGen.PSy`

    :returns: the transformed PSy object.
    :rtype: :py:class:`psyclone.psyGen.PSy`

    '''
    itrans = KernelModuleInlineTrans()

    schedule = psy.invokes.get('invoke_0').schedule

    # Module-Inline all coded kernels in this Schedule
    for kernel in schedule.coded_kernels():
        itrans.apply(kernel)

    return psy
