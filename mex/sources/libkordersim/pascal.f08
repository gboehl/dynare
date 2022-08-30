! Provides a subroutine to build Pascal's triangle
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

module pascal
   implicit none

   type line
      integer, dimension(:), allocatable :: coeffs 
   end type line

   type pascal_triangle
      integer :: d
      type(line), dimension(:), allocatable :: lines     
   end type pascal_triangle

   interface pascal_triangle
      module procedure :: init_pascal_triangle
   end interface pascal_triangle

contains

   ! Fills a pascal_triangle structure with the associated coefficients up to a given order
   type(pascal_triangle) function init_pascal_triangle(d)
      integer, intent(in) :: d
      integer :: i, j
      init_pascal_triangle%d = d
      allocate(init_pascal_triangle%lines(d))
      ! Initializing the first line
      allocate(init_pascal_triangle%lines(1)%coeffs(2))
      init_pascal_triangle%lines(1)%coeffs = 1
      ! Iterating Pascal's triangle formula
      if (d > 1) then
         do i=2,d
            allocate(init_pascal_triangle%lines(i)%coeffs(i+1))
            init_pascal_triangle%lines(i)%coeffs(1) = 1
            init_pascal_triangle%lines(i)%coeffs(i+1) = 1
            do j=2,i
               init_pascal_triangle%lines(i)%coeffs(j) = init_pascal_triangle%lines(i-1)%coeffs(j-1) &
                                                      + init_pascal_triangle%lines(i-1)%coeffs(j)
            end do
         end do
      end if
   end function init_pascal_triangle

   ! Returns ⎛n⎞ stored in pascal_triangle p
   !         ⎝k⎠
   integer function get(k,n,p)
      integer, intent(in) :: k, n
      type(pascal_triangle), intent(in) :: p
      get = p%lines(n)%coeffs(k+1)
   end function get

   ! Returns ⎛        d        ⎞
   !         ⎝ k₁, k₂, ..., kₙ ⎠
   integer function multinomial(k,d,p)
      integer, intent(in) :: k(:), d
      type(pascal_triangle), intent(in) :: p
      integer :: s, i
      s = d 
      multinomial = 1
      i = 1
      do while (s > 0)
         multinomial = multinomial*get(k(i), s, p)         
         s = s-k(i)
         i = i+1
      end do
   end function 

end module pascal

! gfortran -o pascal pascal.f08
! ./pascal
! program test
!    use pascal
!    type(pascal_triangle) :: p
!    integer :: d
!    read *, d
!    p = pascal_triangle(d)
!    do n=1,d 
!       do k=0,n
!          if (k < n) then
!             write (*,'(i2," ")', advance="no") get(k,n,p)
!          else
!             write (*,'(i2," ")') get(k,n,p)
!          end if
!       end do
!    end do

!    d = 3
!    p = pascal_triangle(d)
!    print '(i2)', multinomial([1,2,3], d, p) ! should print 60
!    print '(i2)', multinomial([0,0,0,3], d, p) ! should print 20

! end program test