! Copyright © 2023 Dynare Team
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
!
! Mandatory inputs:
! prhs[1]  Y                  [double] p x T data matrix
! prhs[2]  a                  [double] m x 1 initial mean of the state vector
! prhs[3]  P                  [double] m x m initial variance-covariance matrix
!                                            of the state vector
! prhs[4]  kalman_tol         [double] tolerance parameter (rcond, inversibility
!                                     of the covariance matrix of the prediction
!                                     errors)
! prhs[5]  riccati_tol        [double] tolerance parameter (iteration over the
!                                      Riccati equation). 
! prhs[6]  T                  [double] m x m transition matrix of the state
!                                      equation 
! prhs[7]  Q                  [double] r x r variance-covariance matrix of the
!                                      structural innovations (noise in the 
!                                      state equation).
! prhs[8]  R                  [double] m x r second matrix of the state equation
!                                      relating the structural innovations to
!                                      the state variables.
! prhs[9]  Z                  [double] p x m matrix relating the states to the
!                                      observed variables or vector of indices
!                                      (depending on the value of Zflag).
! Optional inputs:
! prhs[10] Zflag             [integer] equal to 0 if Z is a vector of indices
!                                      targeting the observed variables in the
!                                      state vector (default), equal to 1 if 
!                                      Z is a p x m matrix.
! prhs[11] H                  [double] p x p variance-covariance matrix of the
!                                      measurement errors. If no measurement
!                                      errors set H as a zero scalar (default).
! prhs[12] diffuse_periods   [integer] number of diffuse filter periods in the 
!                                      initialization step (default is 0).   
! prhs[13] presample          [integer] presampling if strictly positive (number
!                                       of initial iterations to be discarded
!                                       when evaluating the likelihood, default
!                                       is 0).
! output:
! plhs[1]  LIK                [double] value of (minus) the likelihood.
! plhs[2]  LIKK               [double] vector containing the density of each
!                                      observation.
subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name="mexFunction")
   use matlab_mex
   use blas
   use lapack
   use iso_c_binding
   implicit none (type, external)

   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs

   type(c_ptr) :: Y_mx, a_mx, P_mx, kalman_tol_mx, &
                 &riccati_tol_mx, presample_mx, T_mx, Q_mx, R_mx, H_mx, &
                 &Z_mx, Zflag_mx, diffuse_periods_mx
   type(c_ptr), dimension(2) :: call_lhs
   type(c_ptr), dimension(11) :: call_rhs
   integer :: presample, diffuse_periods, t, nper, smpl, s
   integer, allocatable, dimension(:) :: indZ
   integer(c_int) :: retval
   integer(blint) :: m, r, p, info, lwork, nb, i, j
   real(real64) :: kalman_tol, riccati_tol, rcond, log_dF, lik
   real(real64), dimension(:,:), contiguous, pointer :: PP, TT, Q, RR, H, Z, Y, K, K_m, iF_m
   real(real64), dimension(:), contiguous, pointer :: a, likk, a_m
   logical :: Zflag, steady_flag, Hflag
   real(real64), allocatable :: v(:), F(:,:), lu(:,:), QQ(:,:),          &
                               &work_rcond(:), work_inv(:), P_next(:,:), &
                               &tmp_m_m(:,:), tmp_m_p(:,:), tmp_m_r(:,:),&
                               &tmp_m_m_prime(:,:), tmp_v(:), a_next(:), &
                               &old_K(:,:), P_iter(:,:), a_iter(:)
   integer(blint), allocatable :: iwork(:), ipiv(:) 

   real(real64), parameter :: pi = 4._real64*DATAN(1._real64)

   ! Check the number of output arguments
   if (nlhs /= 2) then
      call mexErrMsgTxt("Must have 2 outputs")
   end if

   ! Check the consistency and validity of input arguments
   if ((nrhs < 9) .or. (nrhs > 13)) then
      call mexErrMsgTxt("Must have at least 9 inputs and at most 13 inputs")
   end if

   Y_mx = prhs(1)
   a_mx = prhs(2)
   P_mx = prhs(3)
   kalman_tol_mx = prhs(4)
   riccati_tol_mx = prhs(5)
   T_mx = prhs(6)
   Q_mx = prhs(7)
   R_mx = prhs(8)
   Z_mx = prhs(9)

   if (.not. mxIsDouble(Y_mx) .or. mxIsComplex(Y_mx) .or. mxIsSparse(Y_mx)) then
      call mexErrMsgTxt("1st argument (Y) should be a real dense matrix")
   end if
   if (.not. (mxIsDouble(a_mx) .and. ((mxGetM(a_mx) == 1) .or. &
      &(mxGetN(a_mx) == 1))) .or. mxIsComplex(a_mx) .or. mxIsSparse(a_mx)) then
      call mexErrMsgTxt("2nd argument (a) should be a real dense vector")
   end if
   if (.not. mxIsDouble(P_mx) .or. mxIsComplex(P_mx) .or. mxIsSparse(P_mx)) then
      call mexErrMsgTxt("3rd argument (P) should be a real dense matrix")
   end if
   if (.not. (mxIsScalar(kalman_tol_mx) .and. mxIsNumeric(kalman_tol_mx))) then
      call mexErrMsgTxt("4th argument (kalman_tol) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(riccati_tol_mx) .and. mxIsNumeric(riccati_tol_mx))) then
      call mexErrMsgTxt("5th argument (riccati_tol) should be a numeric scalar")
   end if
   if (.not. mxIsDouble(T_mx) .or. mxIsComplex(T_mx) .or. mxIsSparse(T_mx)) then
      call mexErrMsgTxt("6th argument (T) should be a real dense matrix")
   end if
   if (.not. mxIsDouble(Q_mx) .or. mxIsComplex(Q_mx) .or. mxIsSparse(Q_mx)) then
      call mexErrMsgTxt("7th argument (Q) should be a real dense matrix")
   end if
   if (.not. mxIsDouble(R_mx) .or. mxIsComplex(R_mx) .or. mxIsSparse(R_mx)) then
      call mexErrMsgTxt("8th argument (R) should be a real dense matrix")
   end if

   ! Import variables
   a => mxGetPr(a_mx)
   m = size(a, 1, blint)
   p = int(mxGetM(Y_mx), blint)
   nper = int(mxGetN(Y_mx))
   r = int(mxGetM(Q_mx), blint)
   Y(1:p,1:nper) => mxGetPr(Y_mx)
   PP(1:m,1:m) => mxGetPr(P_mx)
   TT(1:m,1:m) => mxGetPr(T_mx)
   Q(1:r,1:r) => mxGetPr(Q_mx)
   RR(1:m,1:r) => mxGetPr(R_mx)
   kalman_tol = mxGetScalar(kalman_tol_mx)
   riccati_tol = mxGetScalar(riccati_tol_mx)

   ! Check the order consistency of input matrices
   if ((mxGetM(P_mx) /= m) .or. (mxGetN(P_mx) /= m) .or.& ! P
      &(mxGetM(T_mx) /= m) .or. (mxGetN(T_mx) /= m) .or.& ! T
      &(mxGetN(Q_mx) /= r) .or.                         & ! Q
      &(mxGetM(R_mx) /= m) .or. (mxGetN(R_mx) /= r)) then ! R
      call mexErrMsgTxt("Input dimension mismatch in (a, Y, P, T, Q, R)")
   end if

   ! Optional inputs
   ! Zflag
   if (nrhs > 9) then
      Zflag_mx = prhs(10)
      if (.not. (mxIsScalar(Zflag_mx) .and. mxIsNumeric(Zflag_mx))) then
         call mexErrMsgTxt("10th argument (Zflag) should be a numeric scalar")
      end if
      Zflag = (mxGetScalar(Zflag_mx) == 1._c_double)
   else
      Zflag = .false.
   end if
   if (Zflag) then
      if (.not. mxIsDouble(Z_mx) .or. mxIsComplex(Z_mx) .or. mxIsSparse(Z_mx)) then
         call mexErrMsgTxt("9th argument (Z) should be a real dense matrix")
      end if
      if ((mxGetM(Z_mx) /= p) .or. (mxGetN(Z_mx) /= m)) then
         call mexErrMsgTxt("Input dimension mismatch in Z")
      end if
      Z(1:p,1:m) => mxGetPr(Z_mx)
   else
      if (.not. (mxIsDouble(Z_mx) .and. ((mxGetM(Z_mx) == 1) .or. &  
         &(mxGetN(Z_mx) == 1))) .or. mxIsComplex(Z_mx) .or.       &
         &mxIsSparse(Z_mx)) then
         call mexErrMsgTxt("9th argument (Z) should be a real dense vector")
      end if
      if ((mxGetM(Z_mx) /= p) .and. (mxGetN(Z_mx) /= p)) then
         call mexErrMsgTxt("Input dimension mismatch in Z")
      end if
      Z(1:p,1:1) => mxGetPr(Z_mx)
      indZ = int(Z(1:p,1))
   end if

   ! H
   if (nrhs > 10) then
      H_mx = prhs(11)
      if ((mxIsScalar(H_mx) .and. mxIsNumeric(H_mx))) then
         Hflag = .false.
      elseif (mxIsDouble(H_mx) .and. .not. (mxIsComplex(H_mx) .or. mxIsSparse(H_mx))) then
         Hflag = .true.
         if ((mxGetM(H_mx) /= p) .or. (mxGetN(H_mx) /= p)) then
            call mexErrMsgTxt("Input dimension mismatch in H")
         end if
         H(1:p,1:p) => mxGetPr(H_mx)
      else
         call mexErrMsgTxt("11th argument (H) should be a real dense matrix or a zero scalar if no measurement error is set")
      end if
   else
      Hflag = .false.
   end if

   ! diffuse_periods
   if (nrhs > 11) then
      diffuse_periods_mx = prhs(12)
      if (.not. (mxIsScalar(diffuse_periods_mx) .and. mxIsNumeric(diffuse_periods_mx))) then
         call mexErrMsgTxt("12th argument (diffuse_periods) should be a numeric scalar")
      end if
      diffuse_periods = int(mxGetScalar(diffuse_periods_mx))
   else
      diffuse_periods = 0
   end if

   ! presample
   if (nrhs > 12) then
      presample_mx = prhs(13)
      if (.not. (mxIsScalar(presample_mx) .and. mxIsNumeric(presample_mx))) then
         call mexErrMsgTxt("13th argument (presample) should be a numeric scalar")
      end if
      presample = int(mxGetScalar(presample_mx))
   else
      presample = 0
   end if

   smpl = nper-diffuse_periods

   ! Density of each observation
   plhs(2) = mxCreateDoubleMatrix(int(smpl, mwSize), 1_mwSize, mxREAL)
   likk(1:smpl) => mxGetPr(plhs(2))

   ! Optimal number of blocks for DGETRI
   nb = ilaenv(1_blint, "DGETRI", " ", p, -1_blint, -1_blint, -1_blint)
   lwork = p*nb

   allocate(v(p), F(p,p), ipiv(p), work_rcond(4*p), work_inv(lwork), iwork(p), &
           &tmp_m_m(m,m), tmp_m_r(m,r), tmp_m_p(m,p), tmp_m_m_prime(m,m), &
           &tmp_v(p), old_K(m,p), K(m,p), a_next(m), QQ(m,m))

   ! Compute RQR'
   ! (i) tmp <- RQ
   call matmul_add("N", "N", 1._real64, RR, Q, 0._real64, tmp_m_r)
   ! (ii) QQ <- tmp*R'
   call matmul_add("N", "T", 1._real64, tmp_m_r, RR, 0._real64, QQ)

   t = diffuse_periods+1
   steady_flag = .false.
   log_dF = huge(log_dF)
   old_K = huge(0._real64)
   a_iter = a
   P_iter = PP
   do 
      if ((t > nper) .or. (steady_flag)) exit
      s = t-diffuse_periods
      ! v <- Y(:,t) - Z*a
      if (Zflag) then
         v = Y(:,t)
         call dgemv("N", p, m, -1._real64, Z, p, a_iter, 1_blint, 1._real64, v, 1_blint)
         ! F <- Z*P*Z' + H
         ! (i) tmp <- P*Z'
         call matmul_add("N", "T", 1._real64, P_iter, Z, 0._real64, tmp_m_p)
         ! (ii) F <- Z*tmp + H
         if (Hflag) then
            F = H
            call matmul_add("N", "N", 1._real64, Z, tmp_m_p, 1._real64, F)
         else
            call matmul_add("N", "N", 1._real64, Z, tmp_m_p, 0._real64, F)
         end if
      else
         v = Y(:,t) - a_iter(indZ)
         if (Hflag) then
            F = P_iter(indZ,indZ) + H
         else
            F = P_iter(indZ,indZ)
         end if
      end if

      ! Compute the reciprocal condition number of F using its LU 
      ! decomposition
      ! (i) LU decomposition of F
      lu = F
      call dgetrf(p, p, lu, p, ipiv, info)
      if (info < 0) then
         call mexErrMsgTxt("Ill-conditioned F!")
      end if
      ! (ii) Reciprocal condition number of F
      call dgecon("1", p, lu, p, norm(F, "1"), rcond, work_rcond, iwork, info)
      if (rcond < kalman_tol) then
         call mexErrMsgTxt("Ill-conditioned F!")
      end if

      ! Compute log(det(F))
      log_dF = 0._real64
      do i=1,p
         log_dF = log_dF+log(abs(lu(i,i)))
      end do
      
      ! Compute the inverse of F using its LU decomposition
      call dgetri(p, lu, p, ipiv, work_inv, lwork, info)

      ! Density of each observation
      call dgemv("N", p, p, 1._real64, lu, p, v, 1_blint, 0._real64, tmp_v, 1_blint)
      likk(s) = log_dF
      do i=1,p
         likk(s) = likk(s)+v(i)*tmp_v(i)
      end do

      ! Compute K
      if (Zflag) then
         ! Compute K = PZ'F^{-1}
         ! (i) tmp <- PZ'
         call matmul_add("N", "T", 1._real64, P_iter, Z, 0._real64, tmp_m_p)
         ! (ii) K <- tmp*F^{-1}
         call matmul_add("N", "N", 1._real64, tmp_m_p, lu, 0._real64, K)
      else
         ! Compute K = P(:,Z)F^{-1} 
         call matmul_add("N", "N", 1._real64, P_iter(:,indZ), lu, 0._real64, K)
      end if

      P_next = QQ

      ! Compute P_{t+1}
      if (Zflag) then
         ! Compute P_{t+1} <- T(I-KZ)PT' + RQR'
         ! (i) Compute tmp' <- I-KZ
         do j=1,m
            do i=1,m
               if (j == i) then
                  tmp_m_m_prime(i,j) = 1._real64
               else
                  tmp_m_m_prime(i,j) = 0._real64
               end if
            end do
         end do
         call matmul_add("N", "N", -1._real64, K, Z, 1._real64, tmp_m_m_prime)
         ! (ii) tmp_m_m <- tmp'*P
         call matmul_add("N", "N", 1._real64, tmp_m_m_prime, P_iter, 0._real64, tmp_m_m)
         ! (iii) tmp_m_m' <- T*tmp_m_m 
         call matmul_add("N", "N", 1._real64, TT, tmp_m_m, 0._real64, tmp_m_m_prime)
         ! (iv) P_next <- tmp_m_m'*T' + P_next
         call matmul_add("N", "T", 1._real64, tmp_m_m_prime, TT, 1._real64, P_next)
      else
         ! Compute P_{t+1} <- T(P-K*P(Z,:))T' + RQR'
         ! (i) tmp_m_m <- P-K*P(Z,:)
         tmp_m_m = P_iter
         call matmul_add("N", "N", -1._real64, K, P_iter(indZ,:), 1._real64, tmp_m_m)
         ! (ii) tmp_m_m' <- T*tmp_m_m
         call matmul_add("N", "N", 1._real64, TT, tmp_m_m, 0._real64, tmp_m_m_prime)
         ! (iii) P_next <- tmp_m_m'*T' + P_next
         call matmul_add("N", "T", 1._real64, tmp_m_m_prime, TT, 1._real64, P_next)
      end if

      ! Compute a_{t+1} = T*(a+Kv)
      ! (i) a <- a+Kv
      call dgemv("N", m, p, 1._real64, K, m, v, 1_blint, 1._real64, a_iter, 1_blint)
      ! (ii) a_next <- T*tmp_v
      call dgemv("N", m, m, 1._real64, TT, m, a_iter, 1_blint, 0._real64, a_next, 1_blint)
      
      ! Check the wedge between gain matrices
      steady_flag = (norm(K-old_K, "M") <= riccati_tol)

      ! Set up next iteration
      old_K = K
      a_iter = a_next
      P_iter = P_next
      t = t+1
   end do

   ! Add observation's densities constants and divide by two.
   likk(1:s) = 0.5_real64*(likk(1:s) + p*log(2*pi))

   ! Call kalman_filter_ss if necessary 
   if (t <= nper) then
      ! call mexPrintf("Calling kalman_filter_ss!")
      call_rhs(1) = Y_mx
      call_rhs(2) = mxCreateDoubleScalar(real(t, c_double))
      call_rhs(3) = mxCreateDoubleScalar(real(nper, c_double))
      call_rhs(4) = mxCreateDoubleMatrix(int(m, mwSize), 1_mwSize, mxREAL)
      a_m => mxGetPr(call_rhs(4))
      a_m = a_iter
      call_rhs(5) = T_mx
      call_rhs(6) = mxCreateDoubleMatrix(int(m, mwSize), int(p, mwSize), mxREAL)
      K_m(1:m,1:p) => mxGetPr(call_rhs(6)) 
      K_m = K
      call_rhs(7) = mxCreateDoubleMatrix(int(p, mwSize), int(p, mwSize), mxREAL)
      iF_m(1:p,1:p) => mxGetPr(call_rhs(7)) 
      iF_m = lu
      call_rhs(8) = mxCreateDoubleScalar(log_dF)
      call_rhs(9) = Z_mx
      call_rhs(10) = mxCreateDoubleScalar(real(p, c_double))
      if (nrhs > 9) then
         call_rhs(11) = Zflag_mx
      else
         call_rhs(11) = mxCreateDoubleScalar(0._c_double)
      end if
      retval = mexCallMATLAB(2_c_int, call_lhs, 11_c_int, call_rhs, "kalman_filter_ss")
      if (retval /= 0_c_int) then
         call mexErrMsgTxt("Error calling kalman_filter_ss!")
      end if
      likk(s+1:smpl) = mxGetPr(call_lhs(2))
   end if

   ! Compute minus the log-likelihood.
   if (presample > diffuse_periods) then
      lik = sum(likk(1+presample-diffuse_periods:smpl))
   else
      lik = sum(likk)
   end if
   plhs(1) = mxCreateDoubleScalar(real(lik, c_double))
      
   ! N.B.: An alternative using Cholesky and the symmetry of matrices goes
   ! (1) P = L*L' (dpotrf)
   ! (2) F <- Z*L*(Z*L)' + H (dsyrk)
   ! F is necessarily symmetric positive-definite in this case
   ! (3) put F in packed format
   ! (4) compute its Cholesky decomposition F = L*L' (dpptrf)
   ! (5) compute the reciprocal condition number via dppcon
   ! (6) Compute P_{t|t} = (I-KZ)P(I-KZ)' + KHK', which ensures symmetry

end subroutine mexFunction