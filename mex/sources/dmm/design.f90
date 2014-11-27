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
#if defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
#if defined(MATLAB_MEX_FILE)
#include "fintrf.h"
#endif

#ifndef MWPOINTER
#define MWPOINTER integer(4)
#endif

#ifndef mwPointer
#define mwPointer MWPOINTER
#endif

SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
  USE MEXINTERFACE
  IMPLICIT NONE

  ! INPUT
  INTEGER(C_INT), INTENT(IN) :: ny,nz,nx,nu,nt
  INTEGER, DIMENSION(6), INTENT(IN) :: ns
  DOUBLE PRECISION, DIMENSION(nt), INTENT(IN) :: theta

  ! OUTPUT
  DOUBLE PRECISION, DIMENSION(ny,max(1,nz),ns(1)), INTENT(OUT) :: c
  DOUBLE PRECISION, DIMENSION(ny,nx,ns(2)), INTENT(OUT) :: H
  DOUBLE PRECISION, DIMENSION(ny,nu,ns(3)), INTENT(OUT) :: G
  DOUBLE PRECISION, DIMENSION(nx,ns(4)), INTENT(OUT) :: a
  DOUBLE PRECISION, DIMENSION(nx,nx,ns(5)), INTENT(OUT) :: F
  DOUBLE PRECISION, DIMENSION(nx,nu,ns(6)), INTENT(OUT) :: R

  ! COMMON
  CHARACTER*200 mfile
  COMMON /M/ mfile

  ! LOCAL
  INTEGER I
  INTEGER(C_INT), DIMENSION(6) :: nsC

  DO I=1,6
     nsC(I) = INT2(ns(I))
  END DO
  CALL designInternal(ny,nz,nx,nu,nsC,nt,theta,TRIM(mfile)//C_NULL_CHAR,c,H,G,a,F,R)

END SUBROUTINE DESIGN
#endif
