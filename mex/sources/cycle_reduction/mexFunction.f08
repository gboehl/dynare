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

! Implements the cyclic reduction algorithm described in
! D.A. Bini, G. Latouche, B. Meini (2002), "Solving matrix polynomial equations arising in queueing problems", Linear Algebra and its Applications 340, pp. 222-244
! D.A. Bini, B. Meini (1996), "On the solution of a nonlinear matrix equation arising in queueing problems", SIAM J. Matrix Anal. Appl. 17, pp. 906-926.
module c_reduction
   use lapack
   use blas
   use matlab_mex
   use ieee_arithmetic
   implicit none (type, external)

contains

   ! Cycle reduction algorithm
   subroutine cycle_reduction(A0, A1, A2, X, cvg_tol, check, info)
      real(real64), dimension(:,:), intent(in) :: A0, A1, A2
      real(real64), dimension(:,:), intent(inout) :: X
      real(real64), intent(in) :: cvg_tol
      logical, intent(in) :: check
      real(c_double), dimension(2), intent(inout) :: info
      
      real(real64), dimension(:,:), allocatable :: A02, &
      Q0, Q2, Ahat, invA1_A02, A1i, A0_tmp, A1_tmp
      integer :: it, n, dn, max_it
      integer(blint) :: info_inv
      integer(blint), dimension(:), allocatable :: ipiv
      real(real64) :: residual, crit
      character(kind=c_char, len=10) :: cvg_tol_str, residual_str

      Ahat = A1
      A1i = A1
      n = size(A0,1)
      dn = 2*n
      allocate(A02(n,dn), ipiv(n), invA1_A02(n,dn), Q0(n,dn), Q2(n,dn))
      A02(:,1:n) = A0
      A02(:,n+1:dn) = A2
      it = 0
      max_it = 300
loop: do
         ! Computing [A0;A2]*(A1\[A0 A2]) 
         A1_tmp = A1i
         invA1_A02 = A02
         call left_divide(A1_tmp, invA1_A02, ipiv, info_inv)
         call matmul_add("N", "N", 1._real64, A02(:,1:n), invA1_A02, 0._real64, Q0)
         call matmul_add("N", "N", 1._real64, A02(:,n+1:dn), invA1_A02, 0._real64, Q2)
         ! Updating A02, A1 and Ahat
         A1i = A1i - Q0(:,n+1:dn) - Q2(:,1:n)
         A02(:,1:n) = -Q0(:,1:n)
         A02(:,n+1:dn) = -Q2(:,n+1:dn)
         Ahat = Ahat-Q2(:,1:n)
         crit = norm(A02(:,1:n),"1")
         ! Checking for stopping conditions
         if (crit < cvg_tol) then
            if (norm(A02(:,n+1:dn), "1") < cvg_tol) then
               exit loop
            end if
         elseif (it == max_it) then
            info(1) = 401._c_double
            info(2) = real(log(norm(A1i,"1")), c_double)
            exit loop
         elseif (ieee_is_nan(crit) .or. (info_inv /= 0_blint)) then
            info(1) = 402._c_double
            info(2) = real(log(norm(A1i,"1")), c_double)
            exit loop
         end if
         it = it + 1
      end do loop

      ! Computing X = -Ahat\A0
      X = -A0
      call left_divide(Ahat, X, ipiv, info_inv)

      ! Checking residuals if necessary
      if (check) then
         ! A1_tmp <- A2*X + A1
         A1_tmp = A1
         call matmul_add("N", "N", 1._real64, A2, X, 1._real64, A1_tmp)
         ! A0_tmp <- A1_tmp*X + A0
         A0_tmp = A0
         call matmul_add("N", "N", 1._real64, A1_tmp, X, 1._real64, A0_tmp)
         residual = norm(A0_tmp, "1")
         if (residual>cvg_tol) then
            info(1) = 403._c_double
            info(2) = log(residual)
            write (cvg_tol_str,"(es8.2)") cvg_tol
            write (residual_str,"(es8.2)") residual
            call mexPrintf("The norm of the residual is "&
                          &// trim(residual_str) // &
                          &", whereas the tolerance criterion is " &
                           // trim(cvg_tol_str) // "." )
         end if
      end if
   end subroutine cycle_reduction

end module c_reduction

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
   use c_reduction
   implicit none (type, external)

   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs

   integer :: i, n
   character(kind=c_char, len=2) :: num2str 
   real(real64) :: cvg_tol
   real(real64), dimension(:,:), pointer,  contiguous :: A0, A1, A2, X
   real(real64), dimension(2) :: info
   logical :: check

   ! 0. Checking the consistency and validity of input arguments
   if (nrhs < 4) then
      call mexErrMsgTxt("Must have at least 4 inputs")
   end if
   if (nrhs > 5) then
      call mexErrMsgTxt("Too many input arguments")
   end if
   if (nlhs > 2) then
      call mexErrMsgTxt("Too many output arguments")
   end if
   
   do i=1,3
      if (.not. (c_associated(prhs(i)) .and. mxIsDouble(prhs(i)) .and. & 
          (.not. mxIsComplex(prhs(i))) .and. (.not. mxIsSparse(prhs(i))))) then
            write (num2str,"(i2)") i
            call mexErrMsgTxt("Argument" // trim(num2str) // " should be a real dense matrix")
      end if
   end do
   if (.not. (c_associated(prhs(4)) .and. mxIsScalar(prhs(4)) .and. &
       mxIsNumeric(prhs(4)))) then
      call mexErrMsgTxt("Argument 4 should be a numeric scalar")
   end if

   check = .false.
   if (nrhs == 5) then
      if (.not. (c_associated(prhs(5)))) then
         call mexErrMsgTxt("Argument 5 should be a Matlab object")
      else
         if (.not. (mxIsEmpty(prhs(5)))) then
            check = .true.
         end if
      end if
   end if

   n = int(mxGetM(prhs(1)))            ! Order of the considered matrices
   if ((n /= mxGetN(prhs(1))) &        ! Number of columns of A0
      &.or. (n /= mxGetM(prhs(2))) &   ! Number of lines of A1
      &.or. (n /= mxGetN(prhs(2))) &   ! Number of columns of A1
      &.or. (n /= mxGetM(prhs(3))) &   ! Number of lines of A2
      &.or. (n /= mxGetN(prhs(3))) &   ! Number of columns of A2
      ) then  
      call mexErrMsgTxt("Input dimension mismatch")
   end if
   
   ! 1. Storing the relevant information in Fortran format
   A0(1:n,1:n) => mxGetPr(prhs(1))
   A1(1:n,1:n) => mxGetPr(prhs(2))
   A2(1:n,1:n) => mxGetPr(prhs(3))
   cvg_tol = mxGetScalar(prhs(4))
   info = [0._c_double,0._c_double]

   plhs(1) = mxCreateDoubleMatrix(int(n, mwSize), int(n, mwSize), mxREAL)
   X(1:n,1:n) => mxGetPr(plhs(1))

   ! 2. Calling the Cycle Reduction algorithm
   call cycle_reduction(A0, A1, A2, X, cvg_tol, check, info)

   ! 3. Editing the information output if necessary
   if (nlhs == 2) then
      if (info(1) == 0.) then
         plhs(2) = mxCreateDoubleScalar(0._c_double)
      else
         plhs(2) = mxCreateDoubleMatrix(1_mwSize, 2_mwSize, mxREAL)
         mxGetPr(plhs(2)) = info
      end if
   end if

end subroutine mexFunction
