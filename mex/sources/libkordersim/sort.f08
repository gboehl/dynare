! Provides a subroutine to sort integer arrays in ascending order
! As the addressed arrays are small, I use the insertion sort algorithm
!
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
! along with Dynare.  If not, see <https://www.gnu.org/licenses/>.

module sort
   implicit none

contains
   subroutine sort_int(l)
      integer, dimension(:), intent(inout) :: l
      integer :: i, j, x
      do i=2,size(l)
         x = l(i)
         j = i
         do while (j > 1 .and. l(j-1) > x)
            l(j) = l(j-1)
            j = j-1
         end do
         l(j) = x
      end do
   end subroutine sort_int

end module sort