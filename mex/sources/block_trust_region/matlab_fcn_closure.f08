! Encapsulates a MATLAB function from ℝⁿ to ℝⁿ
! It must take as 1st argument the point where it is evaluated,
! and return 1 or 2 arguments: value and, optionally, the Jacobian
! The input is copied on entry, the output arguments are copied on exit.
! It may also take extra arguments, which are stored in the module (hence the
! latter is conceptually equivalent to a closure).
!
! Additionally, if x_indices and f_indices are associated, then the matlab_fcn
! procedure only exposes a restricted version of the MATLAB procedure, limited
! to the specified indices specified for x and f. In that case, x_all needs to
! be set in order to give the input values for the indices that are not passed
! to matlab_fcn.

! Copyright © 2019-2020 Dynare Team
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

module matlab_fcn_closure
  use iso_c_binding
  use matlab_mex
  implicit none

  private
  public :: func, extra_args, f_indices, x_indices, x_all, matlab_fcn

  type(c_ptr), pointer :: func => null()
  type(c_ptr), dimension(:), pointer :: extra_args => null()
  integer, dimension(:), pointer :: f_indices => null(), x_indices => null()
  real(real64), dimension(:), pointer :: x_all => null()
contains
  subroutine matlab_fcn(x, fvec, fjac)
    real(real64), dimension(:), intent(in) :: x
    real(real64), dimension(size(x)), intent(out) :: fvec
    real(real64), dimension(size(x), size(x)), intent(out), optional :: fjac

    type(c_ptr), dimension(2) :: call_lhs
    type(c_ptr), dimension(size(extra_args)+2) :: call_rhs
    real(real64), dimension(:), pointer :: x_mat ! Needed to avoid gfortran ICE…
    real(real64), dimension(:), pointer :: fvec_all
    real(real64), dimension(:,:), pointer :: fjac_all
    integer(c_int) :: nlhs
    integer(mwSize) :: n_all

    call_rhs(1) = func
    if (associated(x_all)) then
       n_all = size(x_all)
    else
       n_all = size(x)
    end if
    call_rhs(2) = mxCreateDoubleMatrix(n_all, 1_mwSize, mxREAL)
    x_mat => mxGetPr(call_rhs(2))
    if (associated(x_indices) .and. associated(x_all)) then
       x_mat = x_all
       x_mat(x_indices) = x
    else
       x_mat = x
    end if
    call_rhs(3:) = extra_args

    if (present(fjac)) then
       nlhs = 2
    else
       nlhs = 1
    end if

    ! We use "feval", because it’s the only way of evaluating a function handle through mexCallMATLAB
    if (mexCallMATLAB(nlhs, call_lhs, int(size(call_rhs), c_int), call_rhs, "feval") /= 0) &
         call mexErrMsgTxt("Error calling function to be solved")

    call mxDestroyArray(call_rhs(2))

    fvec_all => mxGetPr(call_lhs(1))
    if (associated(f_indices)) then
       fvec = fvec_all(f_indices)
    else
       fvec = fvec_all
    end if
    call mxDestroyArray(call_lhs(1))

    if (present(fjac)) then
       fjac_all(1:n_all,1:n_all) => mxGetPr(call_lhs(2))
       if (associated(x_indices) .and. associated(f_indices)) then
          fjac = fjac_all(f_indices, x_indices)
       else
          fjac = fjac_all
       end if
       call mxDestroyArray(call_lhs(2))
    end if
  end subroutine matlab_fcn
end module matlab_fcn_closure
