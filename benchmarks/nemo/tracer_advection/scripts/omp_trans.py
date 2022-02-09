# -----------------------------------------------------------------------------
# BSD 3-Clause License
#
# Copyright (c) 2018-2022, Science and Technology Facilities Council
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

'''A simple transformation script for the introduction of OpenMP with PSyclone.

 >>> psyclone -api "nemo" -s ./omp_trans.py tra_adv.F90

This should produce a lot of output, ending with generated Fortran.

'''


def trans(psy):
    ''' Transform a specific Schedule by making all loops
    over levels OpenMP parallel.

    :param psy: the object holding all information on the PSy layer \
                to be modified.
    :type psy: :py:class:`psyclone.psyGen.PSy`

    :returns: the transformed PSy object
    :rtype:  :py:class:`psyclone.psyGen.PSy`

    '''
    from psyclone.psyGen import TransInfo
    from psyclone.nemo import NemoKern
    # Get the transformation we will apply
    ompt = TransInfo().get_trans_name('OMPParallelLoopTrans')
    for invoke in psy.invokes.invoke_list:
        # Get the Schedule of the target routine
        sched = invoke.schedule
        # Apply the OMP transformation to each loop over levels containing
        # a kernel
        for loop in sched.loops():
            kernels = loop.walk(NemoKern)
            if kernels and loop.loop_type == "levels":
                ompt.apply(loop)

    # Return the modified psy object
    return psy
