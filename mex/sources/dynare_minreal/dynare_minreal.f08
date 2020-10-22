! Equivalent to minreal function from MATLAB’s control toolbox
! Essentially a wrapper around SLICOT’s TB01PD.
!
! Takes as first argument a structure with 4 fields:
! — A: n×n real matrix
! − B: n×m real matrix
! − C: p×n real matrix
! − D: p×m real matrix
!
! Returns a similar structure with 4 fields:
! — A: nr×nr real matrix
! − B: nr×m real matrix
! − C: p×nr real matrix
! − D: p×m real matrix
! where nr < n is the size of the minimal state-space system.
!
! An optional second argument sets the tolerance (if omitted, a default value
! is used).

! Copyright © 2021 Dynare Team
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
  use matlab_mex
  use slicot
  implicit none

  type(c_ptr), dimension(*), intent(in), target :: prhs
  type(c_ptr), dimension(*), intent(out) :: plhs
  integer(c_int), intent(in), value :: nlhs, nrhs

  real(real64) :: tol

  integer(c_size_t) :: n, m, p
  integer(slint) :: n_sl, m_sl, p_sl, lda, ldb, ldc, nr, ldwork, info
  type(c_ptr) :: A_mx, B_mx, C_mx, D_mx, Ar_mx, Br_mx, Cr_mx, Dr_mx
  real(real64), dimension(:, :), pointer :: A, B, C, Ar, Br, Cr
  real(real64), dimension(:, :), allocatable :: Atmp, Btmp, Ctmp
  integer(slint), dimension(:), allocatable :: iwork
  real(real64), dimension(:), allocatable :: dwork
  integer(c_int) :: ignore


  !! Process and check arguments

  if (nlhs /= 1 .or. nrhs > 2 .or. nrhs < 1) then
     call mexErrMsgTxt("dynare_minreal: requires 1 or 2 input arguments, and 1 output arguments")
     return
  end if

  if (.not. (mxIsStruct(prhs(1)))) then
     call mexErrMsgTxt("dynare_minreal: first argument (sys) should be a struct")
     return
  end if

  A_mx = mxGetField(prhs(1), 1_mwIndex, "A")
  if (.not. (c_associated(A_mx) .and. mxIsDouble(A_mx) .and. mxGetM(A_mx) == mxGetN(A_mx))) then
     call mexErrMsgTxt("dynare_minreal: first argument (sys) should have a field A containing a square real matrix")
     return
  end if
  n = mxGetM(A_mx)
  A(1:n, 1:n) => mxGetPr(A_mx)

  B_mx = mxGetField(prhs(1), 1_mwIndex, "B")
  if (.not. (c_associated(B_mx) .and. mxIsDouble(B_mx) .and. mxGetM(B_mx) == n)) then
     call mexErrMsgTxt("dynare_minreal: first argument (sys) should have a field B containing a real &
          &matrix with as many lines as sys.A")
     return
  end if
  m = mxGetN(B_mx)
  B(1:n, 1:m) => mxGetPr(B_mx)

  C_mx = mxGetField(prhs(1), 1_mwIndex, "C")
  if (.not. (c_associated(C_mx) .and. mxIsDouble(C_mx) .and. mxGetN(C_mx) == n)) then
     call mexErrMsgTxt("dynare_minreal: first argument (sys) should have a field C containing a real &
          &matrix with as many columns as sys.A")
     return
  end if
  p = mxGetM(C_mx)
  C(1:p, 1:n) => mxGetPr(C_mx)

  D_mx = mxGetField(prhs(1), 1_mwIndex, "D")
  if (.not. (c_associated(D_mx) .and. mxIsDouble(D_mx) .and. mxGetM(D_mx) == p .and. mxGetN(D_mx) == m)) then
     call mexErrMsgTxt("dynare_minreal: first argument (sys) should have a field D containing a real &
          &matrix with as many rows as sys.C and as many columns as sys.B")
     return
  end if

  if (nrhs == 2) then
      if (.not. (mxIsScalar(prhs(2)) .and. mxIsNumeric(prhs(2)))) then
         call mexErrMsgTxt("dynare_minreal: second argument (tol) should be a numeric scalar")
         return
      end if
      tol = mxGetScalar(prhs(2))
   else
      tol = 0._real64 ! Let SLICOT determine the tolerance
   end if


   !! Construct arguments to be passed to the SLICOT function

   n_sl = int(n, slint)
   m_sl = int(m, slint)
   p_sl = int(p, slint)
   lda = max(1_slint, n_sl)
   ldb = max(1_slint, n_sl)
   if (n_sl == 0_slint) then
      ldc = 1_slint
   else
      ldc = max(1_slint, max(m_sl, p_sl))
   end if
   ldwork = max(1_slint, n_sl + max(n_sl, max(3_slint*m_sl, 3_slint*p_sl)))

   allocate(Atmp(n_sl, n_sl), Btmp(ldb, max(m_sl, p_sl)), Ctmp(ldc, n_sl))
   Atmp(1:n_sl, 1:n_sl) = A
   Btmp(1:n_sl, 1:m_sl) = B
   Ctmp(1:p_sl, 1:n_sl) = C

   allocate(iwork(n_sl+max(m_sl,p_sl)), dwork(ldwork))

   call tb01pd("M", "N", n_sl, m_sl, p_sl, A, n_sl, B, m_sl, C, n_sl, nr, tol, iwork, dwork, ldwork, info)


   !! Create output MATLAB structure

   Ar_mx = mxCreateDoubleMatrix(int(nr, c_size_t), int(nr, c_size_t), mxREAL)
   Ar(1:nr, 1:nr) => mxGetPr(Ar_mx)
   Ar = Atmp(1:nr, 1:nr)

   Br_mx = mxCreateDoubleMatrix(int(nr, c_size_t), m, mxREAL)
   Br(1:nr, 1:m) => mxGetPr(Br_mx)
   Br = Btmp(1:nr, 1:m)

   Cr_mx = mxCreateDoubleMatrix(p, int(nr, c_size_t), mxREAL)
   Cr(1:p, 1:nr) => mxGetPr(Cr_mx)
   Cr = Ctmp(1:p, 1:nr)

   Dr_mx = mxDuplicateArray(D_mx)

   plhs(1) = mxCreateStructMatrix(1_mwSize, 1_mwSize)
   ignore = mxAddField(plhs(1), "A")
   call mxSetField(plhs(1), 1_mwIndex, "A", Ar_mx)
   ignore = mxAddField(plhs(1), "B")
   call mxSetField(plhs(1), 1_mwIndex, "B", Br_mx)
   ignore = mxAddField(plhs(1), "C")
   call mxSetField(plhs(1), 1_mwIndex, "C", Cr_mx)
   ignore = mxAddField(plhs(1), "D")
   call mxSetField(plhs(1), 1_mwIndex, "D", Dr_mx)
end subroutine mexFunction
