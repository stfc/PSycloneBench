''' This script estimates the memory usage requirements of the
tracer_advection benchmark. It uses the value of JPI, JPJ and JPK
defined in environment variables as is used by the benchmark, and
outputs the data usage (reads + writes) for a single iteration. This
can be used to determine an estimated memory bandwidth for an execution,
and compare to the cache size(s) of the architecture to tune performance.'''

import os

JPI = int(os.environ.get("JPI"))
JPJ = int(os.environ.get("JPJ"))
JPK = int(os.environ.get("JPK"))

VERBOSE = os.environ.get("VERBOSE", False)
# ELEMENT_SIZE (in bytes) allows for modification of the datatype to be
# other precisions, but assumes double precision by default.
ELEMENT_SIZE = os.environ.get("ELEMENT_SIZE", 8)


def loop_1_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 1 in bytes.
    '''

    size_3d = JPI*JPJ*JPK
    size_2d = JPI*JPJ
    size_1d = JPK

    elements = 3*size_3d + 3*size_2d + 1*size_1d

    byte_used = elements * ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 1 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_2_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 2 in bytes.
    '''

    size_2d = JPI*JPJ

    elements = size_2d

    byte_used = elements * ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 2 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_3_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 3 in bytes.
    '''

    size_2d = JPI*JPJ

    elements = size_2d

    byte_used = elements * ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 3 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_4_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 4 in bytes.
    '''
    size_3d = (JPI-1)*(JPJ-1)*(JPK-1)

    elements = size_3d*5

    byte_used = elements * ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 4 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_5_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 5 in bytes.
    '''
    size_2d = JPI*JPJ

    elements = size_2d

    byte_used = elements * ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 5 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_6_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 6 in bytes.
    '''
    size_2d = JPI*JPJ

    elements = size_2d

    byte_used = elements * ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 6 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_7_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 7 in bytes.
    '''
    size_3d = (JPI-1)*(JPJ-1)*(JPK-1)

    elements = size_3d*4

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 7 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_8_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 8 in bytes.
    '''
    size_3d = (JPI-1)*(JPJ-1)*(JPK-1)

    elements = size_3d*6

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 8 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_9_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 9 in bytes.
    '''
    size_3d = (JPI-2)*(JPJ-2)*(JPK-1)

    elements = size_3d * 8

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 9 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_10_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 10 in bytes.
    '''
    size_3d = (JPI-2)*(JPJ-2)*(JPK-1)

    elements = size_3d * 4

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 10 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_11_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 11 in bytes.
    '''
    size_2d = JPI*JPJ

    elements = size_2d

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 11 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_12_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 12 in bytes.
    '''
    size_2d = JPI*JPJ

    elements = size_2d

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 12 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_13_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 13 in bytes.
    '''
    size_3d = JPI*JPJ*(JPK-2)

    elements = size_3d*3

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 13 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_14_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 14 in bytes.
    '''
    size_2d = JPI*JPJ

    elements = size_2d

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 14 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_15_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 15 in bytes.
    '''
    size_3d = JPI*JPJ*(JPK-2)

    elements = size_3d*2

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 15 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_16_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 16 in bytes.
    '''
    size_3d = JPI*JPJ*(JPK-2)

    elements = size_3d*3

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 16 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_17_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 17 in bytes.
    '''
    size_2d = JPI*JPJ

    elements = size_2d*3

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 17 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_18_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 18 in bytes.
    '''
    size_3d = (JPI-2)*(JPJ-2)*(JPK-1)

    elements = size_3d*5

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 18 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def loop_19_estimate() -> int:
    '''
    :returns: The estimated memory usage of loop 19 in bytes.
    '''
    size_3d = (JPI-2)*(JPJ-2)*(JPK-1)

    elements = size_3d*2

    byte_used = elements*ELEMENT_SIZE

    if VERBOSE:
        print(f"loop 19 uses {byte_used/(1000*1000)} MB.")

    return byte_used


def calculate_memory_usage_per_iteration() -> None:
    ''' Outputs the estimated memory usage (reads and writes) per
    iteration in GB.'''
    byte_used = (
            loop_1_estimate() +
            loop_2_estimate() +
            loop_3_estimate() +
            loop_4_estimate() +
            loop_5_estimate() +
            loop_6_estimate() +
            loop_7_estimate() +
            loop_8_estimate() +
            loop_9_estimate() +
            loop_10_estimate() +
            loop_11_estimate() +
            loop_12_estimate() +
            loop_13_estimate() +
            loop_14_estimate() +
            loop_15_estimate() +
            loop_16_estimate() +
            loop_17_estimate() +
            loop_18_estimate() +
            loop_19_estimate()
    )

    gb_used = byte_used / (1000*1000*1000)

    print(f"Tracer advection model with jpi={JPI}, jpj={JPJ}, and "
          f"jpk={JPK} uses (reads+writes) {gb_used}GB per iteration.")


if __name__ == "__main__":
    calculate_memory_usage_per_iteration()
