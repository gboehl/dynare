! Copyright © 2022 Dynare Team
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

! Implements the logarithmic reduction algorithm described in
! D.A. Bini, G. Latouche, B. Meini (2002), "Solving matrix polynomial equations arising in queueing problems", Linear Algebra and its Applications 340, pp. 222-244

module l_reduction
   use lapack
   use blas
   use matlab_mex
   implicit none

contains

   ! Logarithmic reduction algorithm
   subroutine logarithmic_reduction(A0, A1, A2, G, cvg_tol, check, max_it, info)
      real(real64), dimension(:,:), intent(in) :: A0, A1, A2
      real(real64), dimension(:,:), intent(inout) :: G
      real(real64), intent(in) :: cvg_tol
      logical, intent(in) :: check
      integer, intent(in) :: max_it
      real(c_double), dimension(2), intent(inout) :: info
      
      real(real64), dimension(:,:), allocatable :: tmp_inv, tmp, D02, Id, P, &
     &A1_tmp, A0_tmp, G_new, P_new
      integer :: it, n, dn, i, j
      integer(blint) :: info_inv
      integer(blint), dimension(:), allocatable :: ipiv
      real(real64) :: crit
      character(kind=c_char, len=10) :: cvg_tol_str, residual_str

      n = size(A0,1)
      dn = 2*n
      allocate(D02(n,dn), ipiv(n), Id(n,n), tmp(n,dn), G_new(n,n), P_new(n,n))
      info = [0._c_double,0._c_double]
      ! Set the identity matrix
      do j=1,n
         do i=1,n
            if (i == j) then            
               Id(i,j) = 1.
            else
               Id(i,j) = 0.
            end if
         end do
      end do

      ! Initialization: D02_0
      D02(:,1:n) = A0
      D02(:,n+1:dn) = A2
      tmp_inv = -A1
      call left_divide(tmp_inv, D02, ipiv, info_inv)
      G = D02(:,1:n)
      P = D02(:,n+1:dn)
      it = 0
loop: do
         ! Computing E1_(i+1)
         tmp_inv = Id
         call matmul_add("N", "N", -1._real64, D02(:,1:n), D02(:,n+1:dn), 1._real64, tmp_inv)
         call matmul_add("N", "N", -1._real64, D02(:,n+1:dn), D02(:,1:n), 1._real64, tmp_inv)
         ! Computing matrices [D0^2 D2^2], id est -E0_(i+1) and -E2_(i+1)
         call matmul_add("N", "N", 1._real64, D02(:,1:n), D02(:,1:n), 0._real64, tmp(:,1:n))
         call matmul_add("N", "N", 1._real64, D02(:,n+1:dn), D02(:,n+1:dn), 0._real64, tmp(:,n+1:dn))
         ! Computing D0_(i+1) and D2_(i+1)
         call left_divide(tmp_inv, tmp, ipiv, info_inv)
         ! Computing G_(i+1) = G_i + P_i*D0_(i+1)
         G_new = G
         call matmul_add("N", "N", 1._real64, P, tmp(:,1:n), 1._real64, G_new)
         ! Computing P_(i+1)
         call matmul_add("N", "N", 1._real64, P, tmp(:,n+1:dn), 0._real64, P_new)
         crit = norm(G_new - G,"1")
         ! Checking for stopping conditions
         if (crit < cvg_tol) then
            exit loop
         elseif (it == max_it) then
            info(1) = 411._c_double
            info(2) = real(log(crit), c_double)
            exit loop
         elseif (isnan(crit) .or. (info_inv /= 0_blint)) then
            info(1) = 412._c_double
            info(2) = -1._c_double
            exit loop
         end if
         it = it + 1
         G = G_new
         P = P_new
         D02 = tmp
      end do loop

      ! Checking residuals if necessary
      if (check) then
         ! A1_tmp <- A2*X + A1
         A1_tmp = A1
         call matmul_add("N", "N", 1._real64, A2, G, 1._real64, A1_tmp)
         ! A0_tmp <- A1_tmp*X + A0
         A0_tmp = A0
         call matmul_add("N", "N", 1._real64, A1_tmp, G, 1._real64, A0_tmp)
         crit = norm(A0_tmp, "1")
         if (crit>cvg_tol) then
            info(1) = 413._c_double
            info(2) = real(log(crit), c_double)
            write (cvg_tol_str,"(es8.2)") cvg_tol
            write (residual_str,"(es8.2)") crit
            call mexPrintf("The norm of the residual is "&
                          &// trim(residual_str) // &
                          &", whereas the tolerance criterion is " &
                          // trim(cvg_tol_str) // "." )
         end if
      end if
   end subroutine logarithmic_reduction

end module l_reduction

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
   use l_reduction 
   implicit none   

   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs

   integer :: i, n, max_it
   character(kind=c_char, len=2) :: num2str 
   real(real64) :: cvg_tol
   real(real64), dimension(2) :: info
   real(real64), dimension(:,:), pointer,  contiguous :: A0, A1, A2, X
   logical :: check

   ! 0. Checking the consistency and validity of input arguments
   if (nrhs < 5) then
      call mexErrMsgTxt("Must have at least 5 inputs")
   end if
   if (nrhs > 6) then
      call mexErrMsgTxt("Too many input arguments")
   end if
   if (nlhs > 2) then
      call mexErrMsgTxt("Too many output arguments")
   end if
   
   do i=1,3
      if (.not. (c_associated(prhs(i)) .and. mxIsDouble(prhs(i)) .and. & 
          (.not. mxIsComplex(prhs(i))) .and. (.not. mxIsSparse(prhs(i))))) then
            write (num2str,"(i2)") i
            call mexErrMsgTxt("Argument " // trim(num2str) // " should be a real dense matrix")
      end if
   end do
   do i=4,5
      if (.not. (c_associated(prhs(i)) .and. mxIsScalar(prhs(i)) .and. &
         mxIsNumeric(prhs(i)))) then
            write (num2str,"(i2)") i
            call mexErrMsgTxt("Argument " // trim(num2str) // " should be a numeric scalar")
      end if
   end do

   check = .false.
   if (nrhs == 6) then
      if (.not. (c_associated(prhs(6)))) then
         call mexErrMsgTxt("Argument 6 should be a Matlab object")
      else
         if (.not. (mxIsEmpty(prhs(6)))) then
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
   A2(1:n,1:n) => mxGetPr(prhs(1))
   A1(1:n,1:n) => mxGetPr(prhs(2))
   A0(1:n,1:n) => mxGetPr(prhs(3))
   cvg_tol = mxGetScalar(prhs(4))
   max_it = int(mxGetScalar(prhs(5)))
   info = [0._c_double,0._c_double]

   plhs(1) = mxCreateDoubleMatrix(int(n, mwSize), int(n, mwSize), mxREAL)
   X(1:n,1:n) => mxGetPr(plhs(1))

   ! 2. Calling the Logarithmic Reduction algorithm
   call logarithmic_reduction(A0, A1, A2, X, cvg_tol, check, max_it, info)

   ! 3. Editing the information output if necessary
   if (nlhs == 2) then
      if (info(1) == 0._c_double) then
         plhs(2) = mxCreateDoubleScalar(0._c_double)
      else
         plhs(2) = mxCreateDoubleMatrix(1_mwSize, 2_mwSize, mxREAL)
         mxGetPr(plhs(2)) = info
      end if
   end if

end subroutine mexFunction