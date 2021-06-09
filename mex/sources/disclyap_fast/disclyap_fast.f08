! Solve the discrete Lyapunov Equation (X = G·X·Gᵀ + V) using the Doubling Algorithm
!
! Syntax:
!   [X, error_flag] = disclyap_fast(G, V, tol, check_flag)
!
! Inputs:
!   G             [double]    (n×n) first input matrix
!   V             [double]    (n×n) second input matrix
!   tol           [double]    scalar, tolerance criterion
!   check_flag    [boolean]   if true: check positive-definiteness (optional)
!   max_iter      [integer]   scalar, maximum number of iterations (optional)
!
! Outputs:
!   X             [double]    solution matrix
!   error_flag    [boolean]   true if solution is found, false otherwise (optional)
!
! If check_flag is true, then the code will check if the resulting X
! is positive definite and generate an error message if it is not.
!
! This is a Fortran translation of a code originally written by Joe Pearlman
! and Alejandro Justiniano.

! Copyright © 2020-2021 Dynare Team
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

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
  use iso_fortran_env
  use ieee_arithmetic
  use matlab_mex
  use lapack
  implicit none

  type(c_ptr), dimension(*), intent(in), target :: prhs
  type(c_ptr), dimension(*), intent(out) :: plhs
  integer(c_int), intent(in), value :: nlhs, nrhs

  integer(c_size_t) :: n
  real(real64) :: tol, max_iter
  logical :: check_flag

  real(real64), dimension(:, :), allocatable :: P0, P1, A0, A1, Ptmp
  real(real64) :: matd
  integer :: iter
  integer (blint) :: n_bl

  real(real64), dimension(:, :), pointer :: X

  if (nlhs < 1 .or. nlhs > 2 .or. nrhs < 3 .or. nrhs > 5) then
     call mexErrMsgTxt("disclyap_fast: requires between 3 and 5 input arguments, and 1 or 2 output arguments")
  end if

  n = mxGetM(prhs(1))
  if (.not. mxIsDouble(prhs(1)) .or. mxIsComplex(prhs(1)) &
      .or. .not. mxIsDouble(prhs(2)) .or. mxIsComplex(prhs(2)) &
      .or. mxGetN(prhs(1)) /= n .or. mxGetM(prhs(2)) /= n .or. mxGetN(prhs(2)) /= n) then
     call mexErrMsgTxt("disclyap_fast: first two arguments should be real matrices of the same dimension")
  end if

  if (.not. (mxIsScalar(prhs(3)) .and. mxIsNumeric(prhs(3)))) then
     call mexErrMsgTxt("disclyap_fast: third argument (tol) should be a numeric scalar")
  end if
  tol = mxGetScalar(prhs(3))

  if (nrhs >= 4) then
     if (.not. (mxIsLogicalScalar(prhs(4)))) then
        call mexErrMsgTxt("disclyap_fast: fourth argument (check_flag) should be a logical scalar")
     end if
     check_flag = mxGetScalar(prhs(4)) == 1_c_double
  else
     check_flag = .false.
  end if

  if (nrhs >= 5) then
     if (.not. (mxIsScalar(prhs(5)) .and. mxIsNumeric(prhs(5)))) then
        call mexErrMsgTxt("disclyap_fast: fifth argument (max_iter) should be a numeric scalar")
     end if
     max_iter = int(mxGetScalar(prhs(5)))
  else
     max_iter = 2000
  end if

  ! Allocate and initialize temporary variables
  allocate(P0(n,n), P1(n,n), A0(n,n), A1(n,n), Ptmp(n,n))
  associate (G => mxGetPr(prhs(1)), V => mxGetPr(prhs(2)))
    P0 = reshape(V, [n, n])
    A0 = reshape(G, [n, n])
  end associate
  iter = 1

  n_bl = int(n, blint)
  do
     ! We don't use matmul() for the time being because -fuse-external-blas does
     ! not work as expected under gfortran 8

     ! Ptmp = A0·P0
     call dgemm("N", "N", n_bl, n_bl, n_bl, 1._real64, A0, n_bl, P0, n_bl, 0._real64, Ptmp, n_bl)
     ! P1 = P0+Ptmp·A0ᵀ
     P1 = P0
     call dgemm("N", "T", n_bl, n_bl, n_bl, 1._real64, Ptmp, n_bl, A0, n_bl, 1._real64, P1, n_bl)
     ! A1 = A0·A0
     call dgemm("N", "N", n_bl, n_bl, n_bl, 1._real64, A0, n_bl, A0, n_bl, 0._real64, A1, n_bl)

     matd = maxval(abs(P1-P0))

     P0 = P1
     A0 = A1
     iter = iter + 1

     if (matd <= tol .or. iter == max_iter) exit
  end do

  ! Allocate and set outputs
  plhs(1) = mxCreateDoubleMatrix(n, n, mxREAL)
  X(1:n, 1:n) => mxGetPr(plhs(1))
  if (nlhs > 1) plhs(2) = mxCreateLogicalScalar(.false._mxLogical)

  if (iter == max_iter) then
     X = ieee_value(X, ieee_quiet_nan)
     if (nlhs > 1) then
        call mxDestroyArray(plhs(2))
        plhs(2) = mxCreateLogicalScalar(.true._mxLogical)
     end if
     return
  end if

  X = (P0+transpose(P0))/2._real64

  ! Check that X is positive definite
  if (check_flag .and. nlhs > 1) then
     block
       real(real64), dimension(n, n) :: X2
       integer(blint) :: info
       ! X2=chol(X)
       X2 = X
       call dpotrf("L", n_bl, X2, n_bl, info)
       if (info /= 0) then
          call mxDestroyArray(plhs(2))
          plhs(2) = mxCreateLogicalScalar(.true._mxLogical)
       end if
   end block
  end if
end subroutine mexFunction
