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

module slicot
  use iso_fortran_env

#if defined(MATLAB_MEX_FILE) && __SIZEOF_POINTER__ == 8
  integer, parameter :: slint = int64
  integer, parameter :: sllog = 8 ! Logical kind, gfortran-specific
#else
  integer, parameter :: slint = int32
  integer, parameter :: sllog = 4 ! Logical kind, gfortran-specific
#endif

  interface
     subroutine tb01pd(job, equil, n, m, p, a, lda, b, ldb, c, ldc, &
          nr, tol, iwork, dwork, ldwork, info)
       import :: slint, real64
       character, intent(in) :: job, equil
       integer(slint), intent(in) :: n, m, p, lda, ldb, ldc, ldwork
       real(real64), dimension(*), intent(inout) :: a, b, c
       integer(slint), intent(out) :: nr, info
       real(real64), intent(in) :: tol
       integer(slint), dimension(*), intent(out) :: iwork
       real(real64), dimension(*), intent(out) :: dwork
     end subroutine tb01pd
  end interface
end module slicot
