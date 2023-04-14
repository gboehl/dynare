! Wrapper around MATLAB’s dmperm to compute the Dulmage-Mendelsohn
! decomposition

! Copyright © 2020-2023 Dynare Team
!
! This file is part of Dynare.
!
! Dynare is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Dynare is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

module dulmage_mendelsohn
  use iso_fortran_env
  use matlab_mex
  implicit none (type, external)

  ! Represents a block in the fine DM decomposition
  type :: dm_block
     integer, dimension(:), allocatable :: row_indices, col_indices
  end type dm_block
contains
  subroutine dm_blocks(mat, blocks)
    real(real64), dimension(:, :), intent(in) :: mat
    type(dm_block), dimension(:), allocatable, intent(out) :: blocks

    type(c_ptr), dimension(1) :: call_rhs
    type(c_ptr), dimension(4) :: call_lhs
    real(real64), dimension(:, :), pointer :: mat_mx
    real(real64), dimension(:), pointer :: p, q, r, s
    integer :: i, j

    call_rhs(1) = mxCreateDoubleMatrix(int(size(mat, 1), mwSize), int(size(mat, 2), mwSize), mxREAL)
    mat_mx(1:size(mat,1), 1:size(mat,2)) => mxGetPr(call_rhs(1))
    mat_mx = mat

    if (mexCallMATLAB(4_c_int, call_lhs, 1_c_int, call_rhs, "dmperm") /= 0) &
         call mexErrMsgTxt("Error calling dmperm")

    call mxDestroyArray(call_rhs(1))

    p => mxGetPr(call_lhs(1))
    q => mxGetPr(call_lhs(2))
    r => mxGetPr(call_lhs(3))
    s => mxGetPr(call_lhs(4))

    allocate(blocks(size(r)-1))
    do i = 1, size(r)-1
       allocate(blocks(i)%row_indices(int(r(i+1)-r(i))))
       do j = 1, int(r(i+1)-r(i))
          blocks(i)%row_indices(j) = int(p(j+int(r(i))-1))
       end do
       allocate(blocks(i)%col_indices(int(s(i+1)-s(i))))
       do j = 1, int(s(i+1)-s(i))
          blocks(i)%col_indices(j) = int(q(j+int(s(i))-1))
       end do
    end do

    call mxDestroyArray(call_lhs(1))
    call mxDestroyArray(call_lhs(2))
    call mxDestroyArray(call_lhs(3))
    call mxDestroyArray(call_lhs(4))
  end subroutine dm_blocks
end module dulmage_mendelsohn
