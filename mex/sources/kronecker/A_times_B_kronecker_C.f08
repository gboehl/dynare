! This MEX file computes A·(B⊗C) or A·(B⊗B) without explicitly building B⊗C or
! B⊗B, so that one can consider large matrices B and/or C.

! Copyright © 2007-2023 Dynare Team
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
  use iso_c_binding
  use matlab_mex
  use blas
  implicit none (type, external)

  type(c_ptr), dimension(*), intent(in), target :: prhs
  type(c_ptr), dimension(*), intent(out) :: plhs
  integer(c_int), intent(in), value :: nlhs, nrhs

  integer(c_size_t) :: mA, nA, mB, nB, mC, nC
  real(real64), dimension(:, :), pointer, contiguous :: A, B, C, D

  if (nrhs > 3 .or. nrhs < 2 .or. nlhs /= 1) then
     call mexErrMsgTxt("A_times_B_kronecker_C takes 2 or 3 input arguments and provides 1 output argument")
  end if

  if (.not. mxIsDouble(prhs(1)) .or. mxIsComplex(prhs(1)) .or. mxIsSparse(prhs(1)) &
       .or. .not. mxIsDouble(prhs(2)) .or. mxIsComplex(prhs(2)) .or. mxIsSparse(prhs(2))) then
     call mexErrMsgTxt("A_times_B_kronecker_C: first two arguments should be real dense matrices")
  end if
  mA = mxGetM(prhs(1))
  nA = mxGetN(prhs(1))
  mB = mxGetM(prhs(2))
  nB = mxGetN(prhs(2))
  A(1:mA,1:nA) => mxGetPr(prhs(1))
  B(1:mB,1:nB) => mxGetPr(prhs(2))

  if (nrhs == 3) then
     ! A·(B⊗C) is to be computed.
     if (.not. mxIsDouble(prhs(3)) .or. mxIsComplex(prhs(3)) .or. mxIsSparse(prhs(3))) then
        call mexErrMsgTxt("A_times_B_kronecker_C: third argument should be a real dense matrix")
     end if
     mC = mxGetM(prhs(3))
     nC = mxGetN(prhs(3))
     if (mB*mC /= nA) then
        call mexErrMsgTxt("Input dimension error!")
     end if

     C(1:mC,1:nC) => mxGetPr(prhs(3))

     plhs(1) = mxCreateDoubleMatrix(mA, nB*nC, mxREAL)
     D(1:mA,1:nB*nC) => mxGetPr(plhs(1))

     call full_A_times_kronecker_B_C
  else
     ! A·(B⊗B) is to be computed.
     if (mB*mB /= nA) then
        call mexErrMsgTxt("Input dimension error!")
     end if

     plhs(1) = mxCreateDoubleMatrix(mA, nB*nB, mxREAL)
     D(1:mA,1:nB*nB) => mxGetPr(plhs(1))

     call full_A_times_kronecker_B_B
  end if

contains
  ! Computes D=A·(B⊗C)
  subroutine full_A_times_kronecker_B_C
    integer(c_size_t) :: i, j, ka, kd

    kd = 1
    do j = 1,nB
       ka = 1
       do i = 1,mB
          ! D(:,kd:kd+nC) += B(i,j)·A(:,ka:ka+mC)·C
          call dgemm("N", "N", int(mA, blint), int(nC, blint), int(mC, blint), B(i,j), &
               A(:,ka:ka+mC), int(mA, blint), C, int(mC, blint), 1._real64, &
               D(:,kd:kd+nC), int(mA, blint))
          ka = ka + mC
       end do
       kd = kd + nC
    end do
  end subroutine full_A_times_kronecker_B_C

  ! Computes D=A·(B⊗B)
  subroutine full_A_times_kronecker_B_B
    integer(c_size_t) :: i, j, ka, kd

    kd = 1
    do j = 1,nB
       ka = 1
       do i = 1,mB
          ! D(:,kd:kd+nB) += B(i,j)·A(:,ka:ka+mB)·B
          call dgemm("N", "N", int(mA, blint), int(nB, blint), int(mB, blint), B(i,j), &
               A(:,ka:ka+mB), int(mA, blint), B, int(mB, blint), 1._real64, &
               D(:,kd:kd+nB), int(mA, blint))
          ka = ka + mB
       end do
       kd = kd + nB
    end do
  end subroutine full_A_times_kronecker_B_B
end subroutine mexFunction
