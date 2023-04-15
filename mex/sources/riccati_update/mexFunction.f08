! Copyright © 2022-2023 Dynare Team
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

! Implements Ptmp = T*(P-K*Z*P)*transpose(T)+Q where 
! P is the (r x r) variance-covariance matrix of the state vector
! T is the (r x r) transition matrix of the state vector
! K is the (r x n) gain matrix
! Z is the (n x r) matrix linking observable variables to state variables
! Q is the (r x r) variance-covariance matrix of innovations in the state equation
! and accounting for different properties:
! P is a (symmetric) positive semi-definite matrix
! T can be triangular

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
   use matlab_mex
   use blas
   implicit none (type, external)

   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs

   real(real64), dimension(:,:), pointer, contiguous :: P, T, K, Z, Q, Pnew
   real(real64), dimension(:,:), allocatable :: tmp1, tmp2
   integer :: i, n, r
   character(kind=c_char, len=2) :: num2str 

   ! 0. Checking the consistency and validity of input arguments
   if (nrhs /= 5_c_int) then
      call mexErrMsgTxt("Must have 5 input arguments")
   end if
   if (nlhs > 1_c_int) then
      call mexErrMsgTxt("Too many output arguments")
   end if

   do i=1,5
      if (.not. (c_associated(prhs(i)) .and. mxIsDouble(prhs(i)) .and. & 
          (.not. mxIsComplex(prhs(i))) .and. (.not. mxIsSparse(prhs(i))))) then
            write (num2str,"(i2)") i
            call mexErrMsgTxt("Argument " // trim(num2str) // " should be a real dense matrix")
      end if
   end do

   r = int(mxGetM(prhs(1)))            ! Number of states
   n = int(mxGetN(prhs(3)))            ! Number of observables

   if ((r /= mxGetN(prhs(1))) &        ! Number of columns of P
      &.or. (r /= mxGetM(prhs(2))) &   ! Number of lines of T
      &.or. (r /= mxGetN(prhs(2))) &   ! Number of columns of T
      &.or. (r /= mxGetM(prhs(3))) &   ! Number of lines of K
      &.or. (n /= mxGetM(prhs(4))) &   ! Number of lines of Z
      &.or. (r /= mxGetN(prhs(4))) &   ! Number of columns of Z
      &.or. (r /= mxGetM(prhs(5))) &   ! Number of lines of Q
      &.or. (r /= mxGetN(prhs(5))) &   ! Number of columns of Q
      ) then  
      call mexErrMsgTxt("Input dimension mismatch")
   end if

   ! 1. Storing the relevant information in Fortran format
   P(1:r,1:r) => mxGetPr(prhs(1))
   T(1:r,1:r) => mxGetPr(prhs(2))
   K(1:r,1:n) => mxGetPr(prhs(3))
   Z(1:n,1:r) => mxGetPr(prhs(4))
   Q(1:r,1:r) => mxGetPr(prhs(5))

   plhs(1) = mxCreateDoubleMatrix(int(r, mwSize), int(r, mwSize), mxREAL)
   Pnew(1:r, 1:r) => mxGetPr(plhs(1))

   ! 2. Computing the Riccati update of the P matrix 
   allocate(tmp1(r,r), tmp2(r,r))
   ! Pnew <- Q
   Pnew = Q
   ! tmp1 <- K*Z
   call matmul_add("N", "N", 1._real64, K, Z, 0._real64, tmp1)
   ! tmp2 <- P
   tmp2 = P
   ! tmp2 <- tmp2 - tmp1*P
   call matmul_add("N", "N", -1._real64, tmp1, P, 1._real64, tmp2)
   ! tmp1 <- T*tmp2
   call matmul_add("N", "N", 1._real64, T, tmp2, 0._real64, tmp1)
   ! Pnew <- tmp1*T' + Pnew
   call matmul_add("N", "T", 1._real64, tmp1, T, 1._real64, Pnew)

end subroutine mexFunction
