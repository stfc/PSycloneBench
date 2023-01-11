# -----------------------------------------------------------------------------
# BSD 3-Clause License
#
# Copyright (c) 2018-2021, Science and Technology Facilities Council.
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
directives to NEMO style code.  In order to use
it you must first install PSyclone. See README.md in the top-level
psyclone directory.

Once you have psyclone installed, this may be used by doing:

 $ psyclone -api nemo -s ./kernels_trans.py some_source_file.f90

This should produce a lot of output, ending with generated
Fortran. Note that the Fortran source files provided to PSyclone must
have already been preprocessed (if required).

The transformation script attempts to insert Kernels directives at the
highest possible location(s) in the schedule tree (i.e. to enclose as
much code as possible in each Kernels region). However, due to
limitations in the PGI compiler, we must take care to exclude certain
nodes (such as If blocks) from within Kernel regions. If a proposed
region is found to contain such a node (by the ``valid_kernel``
routine) then the script moves a level down the tree and then repeats
the process of attempting to create the largest possible Kernel
region.

'''

from __future__ import print_function
from psyclone.domain.nemo.transformations import NemoAllArrayRange2LoopTrans
from psyclone.errors import InternalError
from psyclone.psyir.nodes import Assignment, CodeBlock, Call, Literal, Loop, \
    ACCLoopDirective
from psyclone.transformations import ACCLoopTrans, TransformationError, \
    ACCKernelsTrans
from utils import valid_kernel


COLLAPSE_LOOPS = False
EXPAND_IMPLICIT_LOOPS = False

# Get the PSyclone transformations we will use
ARRAY_RANGE_TRANS = NemoAllArrayRange2LoopTrans()


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
        sched.view()

        if EXPAND_IMPLICIT_LOOPS:
            # Transform any array assignments (Fortran ':' notation)
            # into loops.
            for assignment in sched.walk(Assignment):
                ARRAY_RANGE_TRANS.apply(assignment)

        add_kernels(sched.children)
        sched.view()


def add_kernels(children, default_present=True):
    '''
    Walks through the PSyIR inserting OpenACC KERNELS directives at as
    high a level as possible.

    :param children: list of sibling Nodes in PSyIR that are candidates for \
                     inclusion in an ACC KERNELS region.
    :type children: list of :py:class:`psyclone.psyir.nodes.Node`
    :param bool default_present: whether or not to supply the \
                          DEFAULT(PRESENT) clause to ACC KERNELS directives.

    '''
    if not children:
        return

    node_list = []
    for child in children[:]:
        # Can this node be included in a kernels region?
        if not valid_kernel(child):
            try_kernels_trans(node_list, default_present)
            node_list = []
            # It can't so go down a level and try again
            add_kernels(child.children)
        else:
            node_list.append(child)
    try_kernels_trans(node_list, default_present)


def try_kernels_trans(nodes, default_present):
    '''
    Attempt to enclose the supplied list of nodes within a kernels
    region. If the transformation fails then the error message is
    reported but execution continues.

    :param nodes: list of Nodes to enclose within a Kernels region.
    :type nodes: list of :py:class:`psyclone.psyir.nodes.Node`
    :param bool default_present: whether or not to supply the \
                          DEFAULT(PRESENT) clause to ACC KERNELS directives.

    '''
    if not nodes:
        return
    try:
        ACCKernelsTrans().apply(nodes, {"default_present": default_present})
        if COLLAPSE_LOOPS:
            collapse_loops(nodes)
    except (TransformationError, InternalError) as err:
        print(f"Failed to transform nodes: {nodes}")
        print(f"Error was: {err}")


def collapse_loops(nodes):
    '''
    Searches the supplied list of nodes and applies an ACC LOOP COLLAPSE(2)
    directive to any perfectly-nested lon-lat loops.

    :param nodes: list of nodes to search for loops.
    :type nodes: list of :py:class:`psyclone.psyir.nodes.Node`

    '''
    loop_trans = ACCLoopTrans()

    for node in nodes:
        loops = node.walk(Loop)
        for loop in loops:
            if loop.ancestor(ACCLoopDirective):
                # We've already transformed a parent Loop so skip this one.
                continue
            loop_options = {}
            # We put a COLLAPSE(2) clause on any perfectly-nested lon-lat
            # loops that have a Literal value for their step. The latter
            # condition is necessary to avoid compiler errors with 20.7.
            if loop.loop_type == "lat" and \
               isinstance(loop.step_expr, Literal) and \
               isinstance(loop.loop_body[0], Loop) and \
               loop.loop_body[0].loop_type == "lon" and \
               isinstance(loop.loop_body[0].step_expr, Literal) and \
               len(loop.loop_body.children) == 1:
                loop_options["collapse"] = 2
            if loop_options:
                loop_trans.apply(loop, loop_options)
