! Copyright © 2019 Dynare Team
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
! along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
  use iso_fortran_env
  use iso_c_binding
  use dulmage_mendelsohn
  use matlab_mex
  use matlab_fcn_closure
  use trust_region
  implicit none

  type(c_ptr), dimension(*), intent(in), target :: prhs
  type(c_ptr), dimension(*), intent(out) :: plhs
  integer(c_int), intent(in), value :: nlhs, nrhs

  real(real64), dimension(:), allocatable, target :: x
  type(dm_block), dimension(:), allocatable, target :: blocks
  integer :: info, i
  real(real64), parameter :: tolf = 1e-6_real64
  real(real64), dimension(:), allocatable :: fvec
  real(real64), dimension(:,:), allocatable :: fjac
  logical :: debug
  character(len=80) :: debug_msg

  if (nrhs < 3 .or. nlhs /= 2) then
     call mexErrMsgTxt("Must have at least 3 inputs and exactly 2 outputs")
     return
  end if

  if (.not. ((mxIsChar(prhs(1)) .and. mxGetM(prhs(1)) == 1) .or. mxIsClass(prhs(1), "function_handle"))) then
     call mexErrMsgTxt("First argument should be a string or a function handle")
     return
  end if

  if (.not. (mxIsDouble(prhs(2)) .and. (mxGetM(prhs(2)) == 1 .or. mxGetN(prhs(2)) == 1))) then
     call mexErrMsgTxt("Second argument should be a real vector")
     return
  end if

  if (.not. (mxIsLogicalScalar(prhs(3)))) then
     call mexErrMsgTxt("Third argument should be a logical scalar")
     return
  end if

  func => prhs(1)
  debug = mxGetScalar(prhs(3)) == 1._c_double
  extra_args => prhs(4:nrhs)
  associate (x_mat => mxGetPr(prhs(2)))
    allocate(x(size(x_mat)))
    x = x_mat
  end associate

  allocate(fvec(size(x)))
  allocate(fjac(size(x), size(x)))

  ! Compute block decomposition
  nullify(x_indices, f_indices, x_all)
  call matlab_fcn(x, fvec, fjac)
  call dm_blocks(fjac, blocks)

  if (debug) then
     write (debug_msg, "('DYNARE_SOLVE (solve_algo=13|14): number of blocks = ', i0)") size(blocks)
     call mexPrintf_trim_newline(debug_msg)
  end if

  ! Solve the system, starting from bottom-rightmost block
  x_all => x
  do i = size(blocks),1,-1
     if (debug) then
        write (debug_msg, "('DYNARE_SOLVE (solve_algo=13|14): solving block ', i0, ' of size ', i0)") &
             i, size(blocks(i)%col_indices)
        call mexPrintf_trim_newline(debug_msg)
     end if

     block
       real(real64), dimension(size(blocks(i)%col_indices)) :: x_block
       x_indices => blocks(i)%col_indices
       f_indices => blocks(i)%row_indices
       x_block = x(x_indices)
       call trust_region_solve(x_block, matlab_fcn, info, tolf = tolf)
       x(x_indices) = x_block
     end block
  end do

  ! Verify that we have a solution
  ! Note that here we use the ∞-norm, while trust region uses 2-norm; otherwise
  ! this check would almost always fail (because the 2-norm of the full fvec is
  ! larger than the 2-norm of its sub-vectors)
  ! If the check fails, this normally means that the block decomposition was
  ! incorrect (because some element of the Jacobian was numerically zero at the
  ! guess value while not being symbolically zero)
  nullify(x_indices, f_indices, x_all)
  call matlab_fcn(x, fvec)
  if (maxval(abs(fvec)) > tolf) then
     if (debug) &
          call mexPrintf_trim_newline("DYNARE_SOLVE (solve_algo=13|14): residuals still too large, solving for the whole model")
     call trust_region_solve(x, matlab_fcn, info, tolf = tolf)
  else
     info = 1
  end if

  plhs(1) = mxCreateDoubleMatrix(int(size(x, 1), mwSize), 1_mwSize, mxREAL)
  mxGetPr(plhs(1)) = x
  if (info == 1) then
     plhs(2) = mxCreateDoubleScalar(0._c_double)
  else
     plhs(2) = mxCreateDoubleScalar(1._c_double)
  end if
end subroutine mexFunction
