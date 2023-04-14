! Wrapper around LAPACK’s dgges (generalized Schur decomposition) that gives a
! better access to error conditions than does MATLAB’s qz.
!
! Syntax:
!   [ss, tt, zz, sdim, eigval, info] = mjdgges(e, d, qz_criterium, zhreshold)
!
! Inputs:
!   e            [double]  real square (n×n) matrix
!   d            [double]  real square (n×n) matrix
!   qz_criterium [double]  scalar (of the form 1+ε)
!   zhreshold    [double]  used for detecting eigenvalues too close to 0÷0
!
! Outputs:
!   ss           [double]  (n×n) quasi-triangular matrix
!   tt           [double]  (n×n) quasi-triangular matrix
!   zz           [double]  (n×n) orthogonal matrix
!   sdim         [integer] scalar, number of stable eigenvalues
!   eigval       [complex] (n×1) vector of generalized eigenvalues
!   info         [integer] scalar, error code of dgges (or 30 if eigenvalue close to 0÷0)

! Copyright © 2006-2023 Dynare Team
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

#include "defines.F08"

module select_fct_mod
  use iso_fortran_env
  implicit none (type, external)

  real(real64) :: criterium
contains
  logical(bllog) function select_fct(alpha_r, alpha_i, beta)
    use blas

    real(real64), intent(in) :: alpha_r, alpha_i, beta

    select_fct = alpha_r**2 + alpha_i**2 < criterium**2 * beta**2
  end function select_fct
end module select_fct_mod

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
  use iso_fortran_env
  use iso_c_binding
  use select_fct_mod
  use matlab_mex
  use lapack
  implicit none (type, external)

  type(c_ptr), dimension(*), intent(in), target :: prhs
  type(c_ptr), dimension(*), intent(out) :: plhs
  integer(c_int), intent(in), value :: nlhs, nrhs

  integer(c_size_t) :: n
  real(real64) :: zhreshold
  integer(blint) :: n_bl, lwork, info_bl, sdim_bl
  real(real64), dimension(:), allocatable :: alpha_r, alpha_i, beta, work
  logical(bllog), dimension(:), allocatable :: bwork
  ! The pointers used in the LAPACK call are marked as contiguous, to
  ! avoid temporary copies beforehand.
  real(real64), dimension(:), pointer, contiguous :: s, t, z, info, sdim, vsl
#if MX_HAS_INTERLEAVED_COMPLEX
  complex(real64), dimension(:), pointer :: gev
#else
  real(real64), dimension(:), pointer :: gev_r, gev_i
#endif

  if (nrhs < 2 .or. nrhs > 4 .or. nlhs /= 6) then
     call mexErrMsgTxt("MJDGGES: takes 2, 3 or 4 input arguments and exactly 6 output arguments.")
  end if

  n = mxGetM(prhs(1))
  if (.not. mxIsDouble(prhs(1)) .or. mxIsComplex(prhs(1)) .or. mxIsSparse(prhs(1)) &
      .or. .not. mxIsDouble(prhs(2)) .or. mxIsComplex(prhs(2)) .or. mxIsSparse(prhs(2)) &
      .or. mxGetN(prhs(1)) /= n .or. mxGetM(prhs(2)) /= n .or. mxGetN(prhs(2)) /= n) then
     call mexErrMsgTxt("MJDGGES: first two arguments should be real dense matrices of the same dimension")
  end if

  ! Set criterium for stable eigenvalues
  if (nrhs >= 3 .and. mxGetM(prhs(3)) > 0) then
     if (.not. (mxIsScalar(prhs(3)) .and. mxIsNumeric(prhs(3)))) then
        call mexErrMsgTxt("MJDGGES: third argument (qz_criterium) should be a numeric scalar")
     end if
     criterium = mxGetScalar(prhs(3))
  else
     criterium = 1_real64 + 1e-6_real64
  end if

  ! set criterium for 0/0 generalized eigenvalues */
  if (nrhs == 4 .and. mxGetM(prhs(4)) > 0) then
     if (.not. (mxIsScalar(prhs(4)) .and. mxIsNumeric(prhs(4)))) then
        call mexErrMsgTxt("MJDGGES: fourth argument (zhreshold) should be a numeric scalar")
     end if
     zhreshold = mxGetScalar(prhs(4))
  else
     zhreshold = 1e-6_real64
  end if

  plhs(1) = mxCreateDoubleMatrix(n, n, mxREAL)
  plhs(2) = mxCreateDoubleMatrix(n, n, mxREAL)
  plhs(3) = mxCreateDoubleMatrix(n, n, mxREAL)
  plhs(4) = mxCreateDoubleMatrix(1_mwSize, 1_mwSize, mxREAL)
  plhs(5) = mxCreateDoubleMatrix(n, 1_mwSize, mxCOMPLEX)
  plhs(6) = mxCreateDoubleMatrix(1_mwSize, 1_mwSize, mxREAL)

  s => mxGetPr(plhs(1))
  t => mxGetPr(plhs(2))
  sdim => mxGetPr(plhs(4))
#if MX_HAS_INTERLEAVED_COMPLEX
  gev => mxGetComplexDoubles(plhs(5))
#else
  gev_r => mxGetPr(plhs(5))
  gev_i => mxGetPi(plhs(5))
#endif
  info => mxGetPr(plhs(6))
  z => mxGetPr(plhs(3))
  vsl => null()

  ! Copy input matrices, since we can’t modify them
  associate (a => mxGetPr(prhs(1)), b => mxGetPr(prhs(2)))
    s = a
    t = b
  end associate

  n_bl = int(n, blint)
  lwork = 16*n_bl + 16
  allocate(alpha_r(n_bl), alpha_i(n_bl), beta(n_bl), bwork(n_bl), work(lwork))

  call dgges("N", "V", "S", select_fct, n_bl, s, n_bl, t, n_bl, sdim_bl, &
       alpha_r, alpha_i, beta, vsl, n_bl, z, n_bl, work, lwork, bwork, info_bl)

  info = info_bl
  sdim = sdim_bl

#if MX_HAS_INTERLEAVED_COMPLEX
  where (alpha_i == 0_real64 .and. beta == 0_real64)
     gev = alpha_r / beta
  elsewhere
     gev = cmplx(alpha_r, alpha_i, real64) / beta
  end where
#else
  gev_r = alpha_r / beta
  where (alpha_i == 0_real64 .and. beta == 0_real64)
     gev_i = 0_real64
  elsewhere
     gev_i = alpha_i / beta
  end where
#endif

  ! If the ratio of some eigenvalue is too close to 0/0, return specific
  ! error number (only if no other error)
  if (any(abs(alpha_r) <= zhreshold .and. abs(beta) <= zhreshold) .and. info_bl == 0) &
       info = 30
end subroutine mexFunction
