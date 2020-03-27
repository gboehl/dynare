! Copyright © 2006-2020 Dynare Team
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

#include "defines.F08"

module select_fct_mod
  use iso_fortran_env
  implicit none

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
  implicit none

  type(c_ptr), dimension(*), intent(in), target :: prhs
  type(c_ptr), dimension(*), intent(out) :: plhs
  integer(c_int), intent(in), value :: nlhs, nrhs

  integer(c_size_t) :: m1, n1, m2, n2
  real(real64) :: zhreshold
  integer(blint) :: n_bl, lwork, info_bl, sdim_bl
  real(real64), dimension(:), allocatable :: alpha_r, alpha_i, beta, work
  logical(bllog), dimension(:), allocatable :: bwork
  real(real64), dimension(:), pointer :: s, t, z, info, sdim, vsl
#if MX_HAS_INTERLEAVED_COMPLEX
  complex(real64), dimension(:), pointer :: gev
#else
  real(real64), dimension(:), pointer :: gev_r, gev_i
#endif
#if defined(MATLAB_MEX_FILE) && MATLAB_VERSION < 0x0904
  ! Workaround for MKL bug, see below
  real(real64), dimension(:), allocatable, target :: vsl_target
#endif

  if (nrhs < 2 .or. nrhs > 4 .or. nlhs /= 6) then
     call mexErrMsgTxt("MJDGGES: takes 2, 3 or 4 input arguments and exactly 6 output arguments.")
     return
  end if

  m1 = mxGetM(prhs(1))
  n1 = mxGetN(prhs(1))
  m2 = mxGetM(prhs(2))
  n2 = mxGetN(prhs(2))

  if (.not. mxIsDouble(prhs(1)) .or. mxIsComplex(prhs(1)) &
      .or. .not. mxIsDouble(prhs(2)) .or. mxIsComplex(prhs(2)) &
      .or. m1 /= n1 .or. m2 /= n1 .or. m2 /= n2) then
     call mexErrMsgTxt("MJDGGES requires two square real matrices of the same dimension.")
     return
  end if

  ! Set criterium for stable eigenvalues
  if (nrhs >= 3 .and. mxGetM(prhs(3)) > 0) then
     associate (crit_arg => mxGetPr(prhs(3)))
       criterium = crit_arg(1)
     end associate
  else
     criterium = 1_real64 + 1e-6_real64
  end if

  ! set criterium for 0/0 generalized eigenvalues */
  if (nrhs == 4 .and. mxGetM(prhs(4)) > 0) then
     associate (zhresh_arg => mxGetPr(prhs(4)))
       zhreshold = zhresh_arg(1)
     end associate
  else
     zhreshold = 1e-6_real64
  end if

  plhs(1) = mxCreateDoubleMatrix(n1, n1, mxREAL)
  plhs(2) = mxCreateDoubleMatrix(n1, n1, mxREAL)
  plhs(3) = mxCreateDoubleMatrix(n1, n1, mxREAL)
  plhs(4) = mxCreateDoubleMatrix(1_mwSize, 1_mwSize, mxREAL)
  plhs(5) = mxCreateDoubleMatrix(n1, 1_mwSize, mxCOMPLEX)
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
#if defined(MATLAB_MEX_FILE) && MATLAB_VERSION < 0x0904
  ! The left Schur vectors (VSL) are normally not computed, since JOBVSL="N".
  ! But old MKL versions (at least the one shipped with MATLAB R2009b/7.9,
  ! which is MKL 10.1) are buggy, and passing nullptr for VSL leads to a crash.
  ! Hence we need to allocate space for it.
  ! The bug seems to be fixed in MATLAB R2010a/7.10 (MKL 10.2), but we use the
  ! workaround for all versions < R2018a/9.4, since those share the same
  ! ABI and hence the same executables.
  allocate(vsl_target(n1*n1))
  vsl => vsl_target
#else
  vsl => null()
#endif

  ! Copy input matrices, since we can’t modify them
  associate (a => mxGetPr(prhs(1)), b => mxGetPr(prhs(2)))
    s = a
    t = b
  end associate

  n_bl = int(n1, blint)
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
