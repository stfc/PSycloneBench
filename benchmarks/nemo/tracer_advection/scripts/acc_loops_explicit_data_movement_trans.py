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
# Authors: S. Siso, STFC Daresbury Lab

''' PSyclone transformation script to insert OpenACC Parallel Loop directives
to the outermost loop that is parallelisable, including implicit loops.'''

from psyclone.psyGen import TransInfo
from psyclone.psyir.transformations import ACCUpdateTrans
from psyclone.transformations import ACCEnterDataTrans
from utils import insert_explicit_loop_parallelism, normalise_loops


def trans(psy):
    ''' Add OpenACC Parallel Loop directive to all loops, including implicit
    ones to target GPU parallelism.

    :param psy: the PSy object which this script will transform.
    :type psy: :py:class:`psyclone.psyGen.PSy`
    :returns: the transformed PSy object.
    :rtype: :py:class:`psyclone.psyGen.PSy`

    '''
    acc_parallel_trans = TransInfo().get_trans_name('ACCParallelTrans')
    acc_loop_trans = TransInfo().get_trans_name('ACCLoopTrans')

    print("Invokes found:")
    for invoke in psy.invokes.invoke_list:
        print(invoke.name)

        normalise_loops(
                invoke.schedule,
                unwrap_array_ranges=True,
                hoist_expressions=True,
        )

        insert_explicit_loop_parallelism(
                invoke.schedule,
                region_directive_trans=acc_parallel_trans,
                loop_directive_trans=acc_loop_trans,
                collapse=True
        )

        ACCEnterDataTrans().apply(invoke.schedule)
        ACCUpdateTrans().apply(invoke.schedule)

    return psy
