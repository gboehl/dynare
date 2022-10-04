! Copyright © 2019-2021 Dynare Team
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

module blas
  use iso_fortran_env
  implicit none

#if defined(MATLAB_MEX_FILE) && __SIZEOF_POINTER__ == 8
  integer, parameter :: blint = int64
  integer, parameter :: bllog = 8 ! Logical kind, gfortran-specific
#else
  integer, parameter :: blint = int32
  integer, parameter :: bllog = 4 ! Logical kind, gfortran-specific
#endif

  interface
     subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
       import :: blint, real64
       character, intent(in) :: transa, transb
       integer(blint), intent(in) :: m, n, k, lda, ldb, ldc
       real(real64), dimension(*), intent(in) :: a, b
       real(real64), intent(in) :: alpha, beta
       real(real64), dimension(*), intent(inout) :: c
     end subroutine dgemm
  end interface

  interface
     subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
       import :: blint, real64
       character, intent(in) :: trans
       integer(blint), intent(in) :: m, n, lda, incx, incy
       real(real64), dimension(*), intent(in) :: a, x
       real(real64), intent(in) :: alpha, beta
       real(real64), dimension(*), intent(inout) :: y
     end subroutine dgemv
  end interface

contains

  ! Updating a matrix using blas DGEMM
  ! C <- alpha*op(A)*op(B) + beta*C
  subroutine matmul_add(opA, opB, alpha, A, B, beta, C)
     character, intent(in) :: opA, opB
     real(real64), dimension(:,:), intent(in) :: A, B      
     real(real64), dimension(:,:), intent(inout) :: C 
     real(real64), intent(in) :: alpha, beta
     integer(blint) :: m, n, k, lda, ldb
     if (opA == "N") then
        m = int(size(A,1), blint)
        k = int(size(A,2), blint)
        lda = m
     else
        m = int(size(A,2), blint)
        k = int(size(A,1), blint)
        lda = k
     end if 
     if (opB == "N") then
        n = int(size(B,2), blint)        
        ldb = k
     else
        n = int(size(B,1), blint)
        ldb = n
     end if
#ifdef DEBUG
     if ( (opA /= "N") .and. (opA /= "T") .and. (opA /= "C") ) then
        print *, "opA must be either N, T or C"
     end if 
     if ( (opB /= "N") .and. (opB /= "T") .and. (opB /= "C") ) then
        print *, "opB must be either N, T or C"
     end if 
     if (((opA == "N") .and. (opB == "N") .and. (size(A,2) /= size(B,1))) .or.&
        &((opA == "N") .and. (opB /= "N") .and. (size(A,2) /= size(B,2))) .or.&
        &((opA /= "N") .and. (opB == "N") .and. (size(A,1) /= size(B,1))) .or.&
        &((opA /= "N") .and. (opB /= "N") .and. (size(A,1) /= size(B,2))))    &
        then
        print *, "Inconsistent number of columns of op(A) and number of rows &
                 &of op(B)"
     end if
     if (m /= size(C,1)) then
        print *, "Inconsistent number of rows of op(A) and number of rows &
                 &of C"
     end if
#endif
     call dgemm(opA, opB, m, n, k, alpha, A, lda, B, &
               ldb, beta, C, m)
  end subroutine

end module blas

module lapack
  use blas
  implicit none

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

  interface
     function dlange(norm, m, n, a, lda, work)
       import :: blint, real64
       character, intent(in) :: norm
       integer(blint), intent(in) :: m, n, lda
       real(real64), dimension(*), intent(in) :: a
       real(real64), dimension(*), intent(out) :: work
       real(real64) :: dlange
     end function dlange
  end interface

contains

  ! Computing inv(A)*B using lapack DGESV
  ! B <- inv(A)*B
  subroutine left_divide(A, B, ipiv, info)
     real(real64), dimension(:,:), intent(inout) :: A
     real(real64), dimension(:,:), intent(inout) :: B
     integer(blint), dimension(:), intent(inout) :: ipiv
     integer(blint), intent(inout) :: info
     integer(blint) :: m, n
#ifdef DEBUG
     if (size(ipiv) /= size(A,1)) then 
        print *, "Inconsistent number of rows of A and size of ipiv"
     end if
     if (size(A,2) /= size(B,1)) then
        print *, "Inconsistent number of columns of A and number of rows of B"
     end if
#endif
     m = int(size(A,1), blint)
     n = int(size(B,2), blint)
     call dgesv(m, n, A, m, ipiv, B, m, info)
  end subroutine left_divide

  ! Matlab-like norm function calling lapack DLANGE routine
  real(real64) function norm(A, c)
     real(real64), dimension(:,:), intent(in) :: A
     character, intent(in) :: c
     real(real64), dimension(:), allocatable :: work
     integer(blint) :: nrow
     nrow = int(size(A,1),blint)
     norm = dlange(c, nrow, int(size(A,2), blint), &
                  &A, nrow, work)
  end function norm

end module lapack
