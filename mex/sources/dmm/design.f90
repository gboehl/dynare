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
SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
  USE MEXINTERFACE
  IMPLICIT NONE
  ! INPUT
  INTEGER ny,nz,nx,nu,ns(6),nt
  DOUBLE PRECISION theta(nt), nsd(6)

  ! OUTPUT
  DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2))
  DOUBLE PRECISION G(ny,nu,ns(3)),a(nx,ns(4))
  DOUBLE PRECISION F(nx,nx,ns(5)),R(nx,nu,ns(6))

  ! COMMON
  CHARACTER*200 mfile
  COMMON /M/ mfile

  ! LOCAL
  mwPointer INPUT(6), OUTPUT(6)
  INTEGER STATUS, I
  CHARACTER(len=200) :: toprint
  integer*4, PARAMETER :: mxREAL = 0

  ! Matlab mex/mx functions
  INTEGER mexCallMATLAB
  mwPointer mxGetPr, mxCreateDoubleScalar, mxCreateDoubleMatrix
  EXTERNAL mxGetPr, mxCreateDoubleScalar, mxCreateDoubleMatrix, mxCopyPtrToReal8, mexCallMATLAB

  DO I=1,6
     nsd(I) = ns(I)
  END DO

  ! Set input values
  INPUT(1) = mxCreateDoubleScalar(ny*1.d0)
  INPUT(2) = mxCreateDoubleScalar(nz*1.d0)
  INPUT(3) = mxCreateDoubleScalar(nx*1.d0)
  INPUT(4) = mxCreateDoubleScalar(nu*1.d0)
  INPUT(5) = mxCreateDoubleMatrix(1, 6, mxREAL)
  CALL mxCopyReal8ToPtr(nsd, mxGetPr(INPUT(5)), 6)
  INPUT(6) = mxCreateDoubleMatrix(1, nt, mxREAL)
  CALL mxCopyReal8ToPtr(theta, mxGetPr(INPUT(6)), nt)

  ! Call design .m function
  STATUS = mexCallMATLAB(6, OUTPUT, 6, INPUT, mfile)
  IF (STATUS .ne. 0) THEN
     WRITE(toprint, '(A,A,A,I2,A)') '\nCall to ',mfile,' failed with code ',STATUS,'\n'
     CALL mexErrMsgTxt(toprint)
  ENDIF

  ! Copy Matlab output into Fortran
  CALL mxCopyPtrToReal8(mxGetPr(OUTPUT(1)), c, ny*max(1,nz)*ns(1))
  CALL mxCopyPtrToReal8(mxGetPr(OUTPUT(2)), H, ny*nx*ns(2))
  CALL mxCopyPtrToReal8(mxGetPr(OUTPUT(3)), G, ny*nu*ns(3))
  CALL mxCopyPtrToReal8(mxGetPr(OUTPUT(4)), a, nx*ns(4))
  CALL mxCopyPtrToReal8(mxGetPr(OUTPUT(5)), F, nx*nx*ns(5))
  CALL mxCopyPtrToReal8(mxGetPr(OUTPUT(6)), R, nx*nu*ns(6))

  ! Free memory
  DO I = 1,6
     CALL mxDestroyArray(INPUT(I))
     CALL mxDestroyArray(OUTPUT(I))
  END DO
END SUBROUTINE DESIGN
#endif
