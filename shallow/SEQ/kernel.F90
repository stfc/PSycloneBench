!-----------------------------------------------------------------------------
! (c) The copyright relating to this information/data is jointly owned by
! the Crown, Met Office and NERC 2013.
! The contribution of STFC in creating this information/data is acknowledged.
!-----------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! NAME
!   kernel_mod
!
! DESCRIPTION
!   The base type for a kernel. Also contains a component routine that returns a
!   launch object for use by the execution engine/code generator. 
!-------------------------------------------------------------------------------
module kernel_mod
use argument_mod
use lfric
implicit none
private

public CELLS, EDGES, VERTICES, FE, DP, arg
public READ, WRITE, READWRITE, INC
public SUM, MIN, MAX
public CG, DG, R

type :: kernel_type
  private
  logical :: no_op

end type kernel_type

!-------------------------------------------------------------------------------
! Expose public types
!-------------------------------------------------------------------------------

public :: kernel_type

  
end module kernel_mod



