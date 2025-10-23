# -----------------------------------------------------------------------------
# BSD 3-Clause License
#
# Copyright (c) 2022-2023, Science and Technology Facilities Council.
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

''' Utilities file to parallelise Nemo code. '''

import os
from typing import List, Union

from psyclone.errors import InternalError
from psyclone.psyir.nodes import (
    Assignment, Directive, CodeBlock, Call, IfBlock, IntrinsicCall, Loop, Node,
    Reference, Return, Routine, Schedule, StructureReference)
from psyclone.psyir.symbols import DataSymbol
from psyclone.psyir.transformations import (
    ACCKernelsTrans, ArrayAssignment2LoopsTrans, HoistLocalArraysTrans,
    HoistLoopBoundExprTrans,
    HoistTrans, Maxval2LoopTrans, OMPMinimiseSyncTrans, ProfileTrans,
    Reference2ArrayRangeTrans, ScalarisationTrans)
from psyclone.transformations import TransformationError

# If routine names contain these substrings then we do not profile them
PROFILING_IGNORE = []


def normalise_loops(
        schedule,
        hoist_local_arrays: bool = True,
        convert_array_notation: bool = True,
        loopify_array_intrinsics: bool = True,
        convert_range_loops: bool = True,
        scalarise_loops: bool = False,
        hoist_expressions: bool = True,
        ):
    ''' Normalise all loops in the given schedule so that they are in an
    appropriate form for the Parallelisation transformations to analyse
    them.

    :param schedule: the PSyIR Schedule to transform.
    :type schedule: :py:class:`psyclone.psyir.nodes.node`
    :param bool hoist_local_arrays: whether to hoist local arrays.
    :param bool convert_array_notation: whether to convert array notation
        to explicit loops.
    :param bool loopify_array_intrinsics: whether to convert intrinsics that
        operate on arrays to explicit loops (currently only maxval).
    :param bool convert_range_loops: whether to convert ranges to explicit
        loops.
    :param scalarise_loops: whether to attempt to convert arrays to scalars
        where possible, default is False.
    :param hoist_expressions: whether to hoist bounds and loop invariant
        statements out of the loop nest.
    '''
    if hoist_local_arrays:
        # Apply the HoistLocalArraysTrans when possible, it cannot be applied
        # to files with statement functions because it will attempt to put the
        # allocate above it, which is not valid Fortran.
        try:
            HoistLocalArraysTrans().apply(schedule)
        except TransformationError:
            pass

    if convert_array_notation:
        # Make sure all array dimensions are explicit
        for reference in schedule.walk(Reference):
            part_of_the_call = reference.ancestor(Call)
            if part_of_the_call:
                if not part_of_the_call.is_elemental:
                    continue
            if isinstance(reference.symbol, DataSymbol):
                try:
                    Reference2ArrayRangeTrans().apply(reference)
                except TransformationError:
                    pass

    if loopify_array_intrinsics:
        for intr in schedule.walk(IntrinsicCall):
            if intr.intrinsic.name == "MAXVAL":
                try:
                    Maxval2LoopTrans().apply(intr)
                except TransformationError as err:
                    print(err.value)

    if convert_range_loops:
        # Convert all array implicit loops to explicit loops
        explicit_loops = ArrayAssignment2LoopsTrans()
        for assignment in schedule.walk(Assignment):
            if assignment.walk(StructureReference):
                continue  # TODO #2951 Fix issues with structure_refs
            try:
                explicit_loops.apply(assignment)
            except TransformationError:
                pass

    if scalarise_loops:
        # Apply scalarisation to every loop. Execute this in reverse order
        # as sometimes we can scalarise earlier loops if following loops
        # have already been scalarised.
        loops = schedule.walk(Loop)
        loops.reverse()
        scalartrans = ScalarisationTrans()
        for loop in loops:
            scalartrans.apply(loop)

    if hoist_expressions:
        # First hoist all possible expressions
        for loop in schedule.walk(Loop):
            try:
                HoistLoopBoundExprTrans().apply(loop)
            except TransformationError:
                pass

        # Hoist all possible assignments (in reverse order so the inner loop
        # constants are hoisted all the way out if possible)
        for loop in reversed(schedule.walk(Loop)):
            for statement in list(loop.loop_body):
                try:
                    HoistTrans().apply(statement)
                except TransformationError:
                    pass


def insert_explicit_loop_parallelism(
        schedule,
        region_directive_trans=None,
        loop_directive_trans=None,
        collapse: bool = True,
        privatise_arrays: bool = False,
        asynchronous_parallelism: bool = False,
        uniform_intrinsics_only: bool = False,
        enable_reductions: bool = False,
        ):
    ''' For each loop in the schedule that doesn't already have a Directive
    as an ancestor, attempt to insert the given region and loop directives.

    :param schedule: the PSyIR Schedule to transform.
    :type schedule: :py:class:`psyclone.psyir.nodes.node`
    :param region_directive_trans: PSyclone transformation that inserts the
        region directive.
    :type region_directive_trans: \
        :py:class:`psyclone.transformation.Transformation`
    :param loop_directive_trans: PSyclone transformation that inserts the
        loop parallelisation directive.
    :type loop_directive_trans: \
        :py:class:`psyclone.transformation.Transformation`
    :param collapse: whether to attempt to insert the collapse clause to as
        many nested loops as possible.
    :param privatise_arrays: whether to attempt to privatise arrays that cause
        write-write race conditions.
    :param asynchronous_parallelism: whether to attempt to add asynchronocity
    to the parallel sections.
    :param uniform_intrinsics_only: if True it prevent offloading loops
        with non-reproducible device intrinsics.
    :param enable_reductions: whether to enable generation of reduction
        clauses automatically.

    '''
    nemo_v4 = os.environ.get('NEMOV4', False)
    if schedule.name == "ts_wgt":
        return  # TODO #2937 WaW dependency incorrectly considered private
    # Add the parallel directives in each loop
    for loop in schedule.walk(Loop):
        if loop.ancestor(Directive):
            continue  # Skip if an outer loop is already parallelised

        opts = {"collapse": collapse, "privatise_arrays": privatise_arrays,
                "verbose": True, "nowait": asynchronous_parallelism,
                "enable_reductions": enable_reductions}

        if uniform_intrinsics_only:
            opts["device_string"] = "nvfortran-uniform"

        routine_name = loop.ancestor(Routine).name

        if ('dyn_spg' in routine_name and len(loop.walk(Loop)) > 2):
            loop.append_preceding_comment(
                "PSyclone: Loop not parallelised because it is in 'dyn_spg' "
                "and is not the inner loop")
            continue

        try:
            # First check that the region_directive is feasible for this region
            if region_directive_trans:
                # TODO psyclone/#3066 - validate *should* accept a single Node
                # but currently has a bug and doesn't so we have to make a
                # list and pass that.
                region_directive_trans.validate([loop], options=opts)

            # If it is, apply the parallelisation directive
            loop_directive_trans.apply(loop, options=opts)

            # And if successful, the region directive on top.
            if region_directive_trans:
                region_directive_trans.apply(loop.parent.parent, options=opts)
        except TransformationError:
            # This loop cannot be transformed, proceed to next loop.
            # The parallelisation restrictions will be explained with a comment
            # associted to the loop in the generated output.
            continue

    # If we are adding asynchronous parallelism then we now try to minimise
    # the number of barriers.
    if asynchronous_parallelism:
        minsync_trans = OMPMinimiseSyncTrans()
        minsync_trans.apply(schedule)


def add_profiling(children: Union[List[Node], Schedule]):
    '''
    Walks down the PSyIR and inserts the largest possible profiling regions
    in place. Code inside functions or that contains directives is excluded.

    :param children: a Schedule or sibling nodes in the PSyIR to which to
        attempt to add profiling regions.

    '''
    if children and isinstance(children, Schedule):
        # If we are given a Schedule, we look at its children.
        children = children.children

    if not children:
        return

    # We do not want profiling calipers inside functions (such as the
    # PSyclone-generated comparison functions).
    parent_routine = children[0].ancestor(Routine)
    if parent_routine and parent_routine.return_symbol:
        return

    node_list = []
    for child in children[:]:
        # Do we want this node to be included in a profiling region?
        if child.walk((Directive, Return)):
            # It contains a directive or return statement so we put what we
            # have so far inside a profiling region.
            add_profile_region(node_list)
            # A node that is not included in a profiling region marks the
            # end of the current candidate region so reset the list.
            node_list = []
            # Now we go down a level and try again without attempting to put
            # profiling below directives or within Assignments
            if isinstance(child, IfBlock):
                add_profiling(child.if_body)
                add_profiling(child.else_body)
            elif not isinstance(child, (Assignment, Directive)):
                add_profiling(child.children)
        else:
            # We can add this node to our list for the current region
            node_list.append(child)
    add_profile_region(node_list)


def add_profile_region(nodes):
    '''
    Attempt to put the supplied list of nodes within a profiling region.

    :param nodes: list of sibling PSyIR nodes to enclose.
    :type nodes: list of :py:class:`psyclone.psyir.nodes.Node`

    '''
    if nodes:
        # Check whether we should be adding profiling inside this routine
        routine_name = nodes[0].ancestor(Routine).name.lower()
        if any(ignore in routine_name for ignore in PROFILING_IGNORE):
            return
        if len(nodes) == 1:
            if isinstance(nodes[0], CodeBlock) and \
               len(nodes[0].get_ast_nodes) == 1:
                # Don't create profiling regions for CodeBlocks consisting
                # of a single statement
                return
            if isinstance(nodes[0], IfBlock) and \
               "was_single_stmt" in nodes[0].annotations and \
               isinstance(nodes[0].if_body[0], CodeBlock):
                # We also don't put single statements consisting of
                # 'IF(condition) CALL blah()' inside profiling regions
                return
        try:
            ProfileTrans().apply(nodes)
        except TransformationError:
            pass


def valid_kernel(node):
    '''
    Whether the sub-tree that has `node` at its root is eligible to be
    enclosed within an OpenACC KERNELS directive.

    :param node: the node in the PSyIR to check.
    :type node: :py:class:`psyclone.psyir.nodes.Node`

    :returns: True if the sub-tree can be enclosed in a KERNELS region.
    :rtype: bool

    '''
    try:
        ACCKernelsTrans().validate(node, {"disable_loop_check": True})
    except TransformationError:
        return False

    return True


def add_kernels(children: list[Node], default_present: bool = True):
    '''
    Walks through the PSyIR inserting OpenACC KERNELS directives at as
    high a level as possible.

    :param children: list of sibling Nodes in PSyIR that are candidates for
                     inclusion in an ACC KERNELS region.
    :param default_present: whether or not to supply the
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


def try_kernels_trans(nodes: list[Node], default_present: bool):
    '''
    Attempt to enclose the supplied list of nodes within a kernels
    region. If the transformation fails then the error message is
    reported but execution continues.

    :param nodes: list of Nodes to enclose within a Kernels region.
    :param default_present: whether or not to supply the
        DEFAULT(PRESENT) clause to ACC KERNELS directives.

    '''
    if not nodes:
        return
    try:
        ACCKernelsTrans().apply(nodes, {"default_present": default_present})
    except (TransformationError, InternalError) as err:
        print(f"Failed to transform nodes: {nodes}")
        print(f"Error was: {err}")
