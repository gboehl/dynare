!
! Copyright (C) 2014 Dynare Team
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
!

MODULE MEXINTERFACE
  USE ISO_C_BINDING
  IMPLICIT NONE

  INTERFACE

     SUBROUTINE mexPrintf(string) BIND(C, NAME="mexPrintf")
       USE ISO_C_BINDING, ONLY: C_CHAR
       CHARACTER(KIND=C_CHAR) :: string(*)
     END SUBROUTINE mexPrintf

     SUBROUTINE mexErrMsgTxt(string) BIND(C, NAME="mexErrMsgTxt")
       USE ISO_C_BINDING, ONLY: C_CHAR
       CHARACTER(KIND=C_CHAR) :: string(*)
     END SUBROUTINE mexErrMsgTxt

     SUBROUTINE designInternal(ny,nz,nx,nu,ns,nt,theta,mfile,c,H,G,a,F,R) BIND(C, NAME="designInternal")
       USE, INTRINSIC :: ISO_C_BINDING
       IMPLICIT NONE
       INTEGER(C_INT), VALUE, INTENT(IN) :: ny, nz, nx, nu, nt
       INTEGER(C_INT), DIMENSION(*), INTENT(IN) :: ns
       REAL(C_DOUBLE), DIMENSION(*), INTENT(IN) :: theta
       CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: mfile
       REAL(C_DOUBLE), DIMENSION(*), INTENT(OUT) :: c, H, G, a, F, R
     END SUBROUTINE designInternal

  END INTERFACE
END MODULE MEXINTERFACE
