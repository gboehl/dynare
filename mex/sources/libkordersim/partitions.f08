! Provides subroutines to manipulate indexes representing elements of
! a partition for a given integer
! i.e. elements p = (α₁,…,αₘ) where each αᵢ ∈ { 0, ..., n-1 }
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

module partitions
   use pascal
   use sort
   implicit none

   ! index represents the aforementioned (α₁,…,αₘ) objects
   type index
      integer, dimension(:), allocatable :: ind 
   end type index

   interface index
      module procedure :: init_index
   end interface index

   ! a dictionary that matches folded indices with folded offsets
   type dict
      integer :: pr ! pointer to the last added element in indices and offsets
      type(index), dimension(:), allocatable :: indices ! list of folded indices
      integer, dimension(:), allocatable :: offsets ! list of the associated offsets
   end type dict

   interface dict
      module procedure :: init_dict
   end interface dict 

   interface operator(/=)
      module procedure :: diff_indices
   end interface operator(/=)

   ! A type to contain the correspondence unfolded and folded offsets i.e.
   ! folded(i) shall contain the folded offset corresponding to the unfolded offset i
   type uf_matching
      type(integer), dimension(:), allocatable :: folded
   end type

contains

   ! Constructor for the index type
   type(index) function init_index(d, ind)
      integer, intent(in) :: d
      integer, dimension(d), intent(in) :: ind
      allocate(init_index%ind(d))
      init_index%ind = ind
   end function init_index

   ! Comparison for the index type. Returns true if the two indices are different
   type(logical) function diff_indices(i1,i2)
      type(index), intent(in) :: i1, i2 
      if (size(i1%ind) /= size(i2%ind) .or. any(i1%ind /= i2%ind)) then
         diff_indices = .true.
      else
         diff_indices = .false.
      end if
   end function diff_indices

   ! Constructor of the dict type
   type(dict) function init_dict(n, d, p)
      integer, intent(in) :: n, d
      type(pascal_triangle), intent(in) :: p
      integer :: size
      size = get(d, n+d-1, p)
      allocate(init_dict%indices(size), init_dict%offsets(size))
      init_dict%pr = 0
   end function init_dict

   ! Count the number of coordinates similar to the the first one for a given index
   type(integer) function get_prefix_length(idx, d)
      integer, intent(in) :: d
      type(index), intent(in) :: idx
      integer :: i
      i = 1
      if (d>1) then
         do while ((i < d) .and. (idx%ind(i+1) == idx%ind(1)))
            i = i+1
         end do
      end if
      get_prefix_length = i
   end function get_prefix_length

   ! Gets the folded index associated with an unfolded index
   type(index) function u_index_to_f_index(idx, d)
      type(index), intent(in) :: idx
      integer, intent(in) :: d
      u_index_to_f_index = index(d, idx%ind)
      call sort_int(u_index_to_f_index%ind)
   end function u_index_to_f_index 

   ! Converts the offset of an unfolded tensor to the associated unfolded tensor index
   ! Note that the index (α₁,…,αₘ) is such that αᵢ ∈ { 0, ..., n-1 }
   ! and the offset is such that j ∈ {1, ..., nᵈ}
   type(index) function u_offset_to_u_index(j, n, d)
      integer, intent(in) :: j, n, d ! offset, number of variables and dimensions respectively
      integer :: i, tmp, r
      allocate(u_offset_to_u_index%ind(d))
      tmp = j-1 ! We substract 1 as j ∈ {1, ..., n} so that tmp ∈ {0, ..., n-1} and our modular operations work
      do i=d,1,-1
         r = mod(tmp, n)
         u_offset_to_u_index%ind(i) = r
         tmp = (tmp-r)/n
      end do
   end function u_offset_to_u_index

   ! Converts a folded tensor index to the associated folded tensor offset
   ! See the explanation in dynare++/tl/cc/tensor.cc for the function FTensor::getOffsetRecurse
   ! Note that the index (α₁,…,αₘ) is such that αᵢ ∈ { 0, ..., n-1 }
   ! and the offset is such that j ∈ {1, ...,  ⎛n+d-1⎞ }
   !                                           ⎝  d  ⎠
   recursive function f_index_to_f_offset(idx, n, d, p) result(j)
      type(index), intent(in) :: idx   ! folded index
      integer, intent(in) :: n, d      ! number of variables and dimensions
      type(pascal_triangle) :: p       ! Pascal's triangle containing the relevant binomial coefficients
      integer :: j, prefix
      type(index) :: tmp
      if (d == 0) then
         j = 1
      else
         prefix = get_prefix_length(idx,d)
         tmp = index(d-prefix, idx%ind(prefix+1:) - idx%ind(1))
         j = get(d, n+d-1, p) - get(d, n-idx%ind(1)+d-1, p) + f_index_to_f_offset(tmp, n-idx%ind(1), d-prefix, p) 
      end if
   end function f_index_to_f_offset
 
   ! Function that searches a value in an array of a given length 
   type(integer) function find(a, v, l)
      integer, intent(in) :: l ! length of the array
      type(index), dimension(l), intent(in) :: a ! array of indices
      type(index) :: v ! element to be found
      integer :: i
      if (l == 0) then
         find = 0
      else
         i = 1
         do while (i <= l .and. a(i) /= v)
            i = i+1
         end do
         if (i == l+1) then
            find = 0
         else
            find = i
         end if
      end if
   end function find

   ! Fills the folded offset array:
   ! folded(i) shall contain the folded offset corresponding to the unfolded offset i
   ! For each unfolded tensor offset
   ! 	(a) compute the associated unfolded index (u_offset_to_u_index)
   ! 	(b) compute the associated folded index (u_index_to_f_index)
   ! 	(c) has the folded offset already been computed ?
   ! 		(i)   If yes, get the corresponding offset
   ! 		(ii)  If no, compute it (f_index_to_f_offset) and store it for reuse and as a result
   subroutine fill_folded_indices(folded, n, d, p) 
      integer, intent(in) :: n, d
      integer, dimension(n**d), intent(inout) :: folded 
      type(pascal_triangle), intent(in) :: p
      type(dict) :: c
      type(index) :: tmp
      integer :: j, found
      c = dict(n, d, p)
      do j=1,n**d
         tmp = u_offset_to_u_index(j,n,d)
         tmp = u_index_to_f_index(tmp, d)
         found = find(c%indices, tmp, c%pr)
         if (found == 0) then
           c%pr = c%pr+1
           c%indices(c%pr) = tmp
           c%offsets(c%pr) = f_index_to_f_offset(tmp,n,d,p)
           folded(j) = c%offsets(c%pr)  
         else
           folded(j) = c%offsets(found) 
         end if
      end do
   end subroutine fill_folded_indices

   end module partitions

! gfortran -o partitions partitions.f08 pascal.f08 sort.f08
! ./partitions
! program test
!    use partitions
!    use pascal
!    implicit none
!    type(index) :: uidx, fidx, i1, i2
!    integer, dimension(:), allocatable :: folded
!    integer :: i, uj, n, d
!    type(pascal_triangle) :: p
!    ! Unfolded indices and offsets
!    ! 0,0,0  1    1,0,0  10   2,0,0  19
!    ! 0,0,1  2    1,0,1  11   2,0,1  20
!    ! 0,0,2  3    1,0,2  12   2,0,2  21
!    ! 0,1,0  4    1,1,0  13   2,1,0  22
!    ! 0,1,1  5    1,1,1  14   2,1,1  23
!    ! 0,1,2  6    1,1,2  15   2,1,2  24
!    ! 0,2,0  7    1,2,0  16   2,2,0  25
!    ! 0,2,1  8    1,2,1  17   2,2,1  26
!    ! 0,2,2  9    1,2,2  18   2,2,2  27  

!    ! Folded indices and offsets
!    ! 0,0,0  1     1,1,1   7      2,2,2   10
!    ! 0,0,1  2     1,1,2   8
!    ! 0,0,2  3     1,2,2   9
!    ! 0,1,1  4
!    ! 0,1,2  5
!    ! 0,2,2  6
 
!    n = 3
!    d = 3
!    uj = 8
!    p = pascal_triangle(n+d-1)   

!    ! u_offset_to_u_index
!    uidx = u_offset_to_u_index(uj,n,d)
!    print '(3i2)', (uidx%ind(i), i=1,d) ! should display 0 2 1

!    ! f_index_to_f_offset
!    fidx = u_index_to_f_index(uidx, d)
!    print '(i2)', f_index_to_f_offset(fidx, n, d, p) ! should display 5

!    ! /=
!    i1 = index(5, (/1,2,3,4,5/))
!    i2 = index(5, (/1,2,3,4,6/))
!    if (i1 /= i2) then
!       print *, "Same!"
!    else
!       print *, "Different!"
!    end if

!    ! fill_folded_indices
!    allocate(folded(n**d))
!    call fill_folded_indices(folded,n,d)
!    print *, "Matching offsets unfolded -> folded"
!    print '(1000i4)', (i, i=1,n**d)
!    print '(1000i4)', (folded(i), i=1,n**d)

! end program test