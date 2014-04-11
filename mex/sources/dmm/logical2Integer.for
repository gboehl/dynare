C
C Copyright (C) 2014 Dynare Team
C
C This file is part of Dynare.
C
C Dynare is free software: you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation, either version 3 of the License, or
C (at your option) any later version.
C
C Dynare is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
C
      FUNCTION logical2Integer(v, n)
      INTEGER, INTENT(IN) :: n
      LOGICAL, DIMENSION(n), INTENT(IN) :: v
      INTEGER, DIMENSION(n), INTENT(OUT) :: logicalToInteger
      DO I=1,n
         IF (v(I)) THEN
            logicalToInteger(I) = 1
         ELSE
            logicalToInteger(I) = 0
         END IF
      END DO
      END FUNCTION
