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
use global_parameters_mod
implicit none
private

public CELLS, EDGES, VERTICES, FE, DP, arg
public READ, WRITE, READWRITE, INC
public SUM, MIN, MAX

!> These quantities should be defined somewhere in the lfric
!! infrastructure but at the moment they are not!
!! \todo Work out where POINTWISE and DOFS should be declared.
integer, public, parameter :: POINTWISE = 2, DOFS = 5

!> The points in the domain that a kernel will update
integer, public, parameter :: INTERNAL_PTS = 0, &
                              EXTERNAL_PTS = 1, &
                              ALL_PTS = 2
!> The type of grid that a kernel is written to operate on.
!! NEMO uses an orthogonal, curvilinear mesh while
!! shallow has a regular mesh (constant spatial resolution).
!! \todo These should probably be declared somewhere else!
integer, public, parameter :: ORTHOGONAL_REGULAR     = 7
integer, public, parameter :: ORTHOGONAL_CURVILINEAR = 8

type :: kernel_type
  private
  logical :: no_op

end type kernel_type

!-------------------------------------------------------------------------------
! Expose public types
!-------------------------------------------------------------------------------

public :: kernel_type

  
end module kernel_mod



