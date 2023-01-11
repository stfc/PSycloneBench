# -----------------------------------------------------------------------------
# BSD 3-Clause License
#
# Copyright (c) 2018-2023, Science and Technology Facilities Council.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# -----------------------------------------------------------------------------
# Authors: R. W. Ford, A. R. Porter and S. Siso, STFC Daresbury Lab

'''A transformation script that seeks to apply OpenACC KERNELS, LOOP and
explicit DATA directives to NEMO style code.  In order to use it you must first
install PSyclone. See README.md in the top-level psyclone directory.

Once you have psyclone installed, this may be used by doing:

 $ psyclone -api nemo -s ./acc_mixed_explicit_data_movement_trans.py <file>

This should produce a lot of output, ending with generated Fortran. Note
that the Fortran source files provided to PSyclone must have already been
preprocessed (if required).

'''

from psyclone.psyir.transformations import ACCUpdateTrans
from psyclone.transformations import ACCEnterDataTrans, ACCLoopTrans
from utils import add_kernels, normalise_loops, \
    insert_explicit_loop_parallelism


def trans(psy):
    '''A PSyclone-script compliant transformation function. Applies
    OpenACC 'kernels', 'loop' and explicit 'data' directives to NEMO code.

    :param psy: The PSy layer object to apply transformations to.
    :type psy: :py:class:`psyclone.psyGen.PSy`
    '''

    print("Invokes found:")
    print("\n".join([str(name) for name in psy.invokes.names]))

    for invoke in psy.invokes.invoke_list:

        sched = invoke.schedule
        if not sched:
            print("Invoke {invoke.name} has no Schedule! Skipping...")
            continue

        # Convert array and range syntax to explicit loops
        normalise_loops(
            invoke.schedule,
            unwrap_array_ranges=True,
            hoist_expressions=True,
        )

        # Add OpenACC Loop directives
        insert_explicit_loop_parallelism(
            invoke.schedule,
            region_directive_trans=None,
            loop_directive_trans=ACCLoopTrans(),
            collapse=True
        )

        # Add OpenACC Kernel directives
        add_kernels(sched.children)

        # Add OpenACC data directives
        ACCEnterDataTrans().apply(invoke.schedule)
        ACCUpdateTrans().apply(invoke.schedule)
