# -----------------------------------------------------------------------------
# BSD 3-Clause License
#
# Copyright (c) 2019-2020, Science and Technology Facilities Council.
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
# Author: A. R. Porter, STFC Daresbury Lab

''' Module containing various utilities to aid in the application of
    OpenACC KERNELS directives to NEMO source. Mainly required to
    workaround the vagaries of the PGI compiler's support for OpenACC.
'''

from __future__ import print_function

from psyclone.errors import InternalError
from psyclone.psyir.nodes import CodeBlock, Call, Literal, Loop, \
    ACCLoopDirective
from psyclone.transformations import ACCLoopTrans, TransformationError, \
    ACCKernelsTrans


def valid_kernel(node):
    '''
    Whether the sub-tree that has `node` at its root is eligible to be
    enclosed within an OpenACC KERNELS directive.

    :param node: the node in the PSyIR to check.
    :type node: :py:class:`psyclone.psyir.nodes.Node`

    :returns: True if the sub-tree can be enclosed in a KERNELS region.
    :rtype: bool

    '''
    excluded_node_types = (CodeBlock, Call)
    if node.walk(excluded_node_types):
        return False
    return True


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
        _, _ = ACCKernelsTrans().apply(nodes,
                                       {"default_present": default_present})
        # collapse_loops(nodes)
    except (TransformationError, InternalError) as err:
        print("Failed to transform nodes: {0}", nodes)
        print("Error was: {0}".format(str(err)))


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
