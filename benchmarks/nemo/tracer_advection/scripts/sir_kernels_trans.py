'''
'''

from __future__ import print_function
from kernels_trans import add_kernels
from sir_trans import make_sir_compliant


def trans(psy):
    '''
    '''
    for invoke in psy.invokes.invoke_list:

        sched = invoke.schedule
        if not sched:
            print("Invoke {0} has no Schedule! Skipping...".
                  format(invoke.name))
            continue

        make_sir_compliant(sched)
        add_kernels(sched.children)
        sched.view()
