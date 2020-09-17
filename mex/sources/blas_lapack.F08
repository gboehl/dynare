! Copyright © 2019-2020 Dynare Team
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

module blas
  use iso_fortran_env

#if defined(MATLAB_MEX_FILE) && __SIZEOF_POINTER__ == 8
  integer, parameter :: blint = int64
  integer, parameter :: bllog = 8 ! Logical kind, gfortran-specific
#else
  integer, parameter :: blint = int32
  integer, parameter :: bllog = 4 ! Logical kind, gfortran-specific
#endif

  interface
     subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
       import :: blint, real64
       character :: trans
       integer(blint), intent(in) :: m, n, lda, incx, incy
       real(real64), dimension(*), intent(in) :: a, x
       real(real64), intent(in) :: alpha, beta
       real(real64), dimension(*), intent(inout) :: y
     end subroutine dgemv
  end interface
end module blas

module lapack
  use blas

  interface
     subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
       import :: blint, real64
       integer(blint), intent(in) :: n, nrhs, lda, ldb
       real(real64), dimension(*), intent(inout) :: a, b
       integer(blint), dimension(*), intent(out) :: ipiv
       integer(blint), intent(out) :: info
     end subroutine dgesv
  end interface

  interface
     subroutine dgges(jobvsl, jobvsr, sort, selctg, n, a, lda, b, ldb, sdim, &
          alphar, alphai, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, bwork, &
          info)
       import :: blint, bllog, real64
       character, intent(in) :: jobvsl, jobvsr, sort
       interface
          logical(bllog) function selctg(alphar, alphai, beta)
            import :: bllog, real64
            real(real64), intent(in) :: alphar, alphai, beta
          end function selctg
       end interface
       integer(blint), intent(in) :: n, lda, ldb, ldvsl, ldvsr, lwork
       real(real64), dimension(*), intent(inout) :: a, b
       real(real64), dimension(*), intent(out) :: alphar, alphai, beta, vsl, vsr, work
       logical(bllog), dimension(*), intent(out) :: bwork
       integer(blint), intent(out) :: sdim, info
     end subroutine dgges
  end interface

  interface
     subroutine dpotrf(uplo, n, a, lda, info)
       import :: blint, real64
       character, intent(in) :: uplo
       integer(blint), intent(in) :: n, lda
       real(real64), dimension(*), intent(inout) :: a
       integer(blint), intent(out) :: info
     end subroutine dpotrf
  end interface
end module lapack
