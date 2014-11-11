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

     INTEGER(4) FUNCTION mexPrintf(toprint) BIND(C, NAME="mexPrintf")
       USE ISO_C_BINDING
       IMPLICIT NONE
       CHARACTER(KIND=C_CHAR), INTENT(IN) :: toprint(*)
     END FUNCTION mexPrintf

     SUBROUTINE mexErrMsgTxt(toprint) BIND(C, NAME="mexErrMsgTxt")
       USE ISO_C_BINDING
       IMPLICIT NONE
       CHARACTER(KIND=C_CHAR), INTENT(IN) :: toprint(*)
     END SUBROUTINE mexErrMsgTxt

  END INTERFACE
END MODULE MEXINTERFACE
