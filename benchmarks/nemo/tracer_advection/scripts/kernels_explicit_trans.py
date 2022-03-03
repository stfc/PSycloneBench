# -----------------------------------------------------------------------------
# BSD 3-Clause License
#
# Copyright (c) 2022, Science and Technology Facilities Council.
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

'''A transformation script that seeks to apply OpenACC KERNELS
and data-movement directives to NEMO-style code.  In order to use
it you must first install PSyclone. See README.md in the top-level
psyclone directory.

Once you have psyclone installed, this may be used by doing:

 $ psyclone -api nemo -s ./kernels_explicit_trans.py some_source_file.f90

This should produce a lot of output, ending with generated
Fortran. Note that the Fortran source files provided to PSyclone must
have already been preprocessed (if required).

Most of the functionality of this transformation script is imported
from kernels_trans.py so see that file for details. This script extends
the one in the latter file by adding directives to control data movement
between CPU and GPU.

'''

from kernels_trans import add_kernels
from psyclone.domain.nemo.transformations import NemoAllArrayRange2LoopTrans
from psyclone.psyir.nodes import Assignment
from psyclone.psyir.transformations import (ACCUpdateTrans,
                                            HoistLocalArraysTrans)
from psyclone.transformations import ACCEnterDataTrans


EXPAND_IMPLICIT_LOOPS = False

# Get the PSyclone transformations we will use
ARRAY_RANGE_TRANS = NemoAllArrayRange2LoopTrans()
EDATA_TRANS = ACCEnterDataTrans()
UPDATE_TRANS = ACCUpdateTrans()
HOIST_TRANS = HoistLocalArraysTrans()


def trans(psy):
    '''A PSyclone-script compliant transformation function. Applies
    OpenACC 'kernels' and 'data' directives to NEMO code.

    :param psy: The PSy layer object to apply transformations to.
    :type psy: :py:class:`psyclone.psyGen.PSy`
    '''

    print("Invokes found:\n{0}\n".format(
        "\n".join([str(name) for name in psy.invokes.names])))

    for invoke in psy.invokes.invoke_list:

        sched = invoke.schedule
        if not sched:
            print("Invoke {0} has no Schedule! Skipping...".
                  format(invoke.name))
            continue

        # Ensure any local, automatic arrays are promoted to module
        # scope so that they can be kept on the GPU.
        if not sched.is_program:
            HOIST_TRANS.apply(sched)

        if EXPAND_IMPLICIT_LOOPS:
            # Transform any array assignments (Fortran ':' notation)
            # into loops.
            for assignment in sched.walk(Assignment):
                ARRAY_RANGE_TRANS.apply(assignment)

        add_kernels(sched.children)
        EDATA_TRANS.apply(sched)

        UPDATE_TRANS.apply(sched, options={"allow-codeblocks": True})

        sched.view()
