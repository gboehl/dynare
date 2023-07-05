! Copyright © 2021-2023 Dynare Team
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
   use iso_fortran_env
   implicit none (type, external)

   ! index represents the aforementioned (α₁,…,αₘ) objects
   type index
      integer, dimension(:), allocatable :: coor 
   end type index

   interface index
      module procedure :: init_index, init_index_vec, init_index_int
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
   

   ! A type to contain the integer partitions up to integer n
   ! with up to n parts
   type partition_triangle
      type(index), dimension(:), allocatable :: partition ! integer partitions
      integer, dimension(:), allocatable :: count ! numbers of equivalent permuted partitions
   end type partition_triangle

contains

   ! Constructors for the index type
   ! Simply allocates the index with the size provided as input
   type(index) function init_index(d)
      integer, intent(in) :: d
      allocate(init_index%coor(d))
   end function init_index

   ! Creates an index with the vector provided as inputs
   type(index) function init_index_vec(ind)
      integer, dimension(:), intent(in) :: ind
      allocate(init_index_vec%coor(size(ind)))
      init_index_vec%coor = ind
   end function init_index_vec

   ! Creates the index with a given size
   ! and fills it with a given integer
   type(index) function init_index_int(d, m)
      integer, intent(in) :: d, m
      integer :: i
      allocate(init_index_int%coor(d))
      do i=1,d
         init_index_int%coor(i) = m
      end do
   end function init_index_int

   ! Operators for the index type
   ! Comparison for the index type. Returns true if the two indices are different
   type(logical) function diff_indices(i1,i2)
      type(index), intent(in) :: i1, i2 
      if (size(i1%coor) /= size(i2%coor) .or. any(i1%coor /= i2%coor)) then
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
         do while ((i < d) .and. (idx%coor(i+1) == idx%coor(1)))
            i = i+1
         end do
      end if
      get_prefix_length = i
   end function get_prefix_length

   ! Gets the folded index associated with an unfolded index
   type(index) function u_index_to_f_index(idx)
      type(index), intent(in) :: idx
      u_index_to_f_index = index(idx%coor)
      call sort_int(u_index_to_f_index%coor)
   end function u_index_to_f_index 

   ! Converts the offset of an unfolded tensor to the associated unfolded tensor index
   ! Note that the index (α₁,…,αₘ) is such that αᵢ ∈ { 0, ..., n-1 }
   ! and the offset is such that j ∈ {1, ..., nᵈ}
   type(index) function u_offset_to_u_index(j, n, d)
      integer, intent(in) :: j, n, d ! offset, number of variables and dimensions respectively
      integer :: i, tmp, r
      allocate(u_offset_to_u_index%coor(d))
      tmp = j-1 ! We substract 1 as j ∈ {1, ..., n} so that tmp ∈ {0, ..., n-1} and our modular operations work
      do i=d,1,-1
         r = mod(tmp, n)
         u_offset_to_u_index%coor(i) = r
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
         tmp = index(idx%coor(prefix+1:) - idx%coor(1))
         j = get(d, n+d-1, p) - get(d, n-idx%coor(1)+d-1, p) + f_index_to_f_offset(tmp, n-idx%coor(1), d-prefix, p) 
      end if
   end function f_index_to_f_offset

   ! Returns the unfolded tensor offset associated with an unfolded tensor index
   ! Written in a recursive way, the unfolded offset off(α₁,…,αₘ) associated with the
   ! index (α₁,…,αₘ) with αᵢ ∈ {1, ..., n} verifies
   ! off(α₁,…,αₘ) = n*off(α₁,…,αₘ₋₁) + αₘ
   integer function u_index_to_u_offset(idx, n, d)
      type(index), intent(in) :: idx   ! unfolded index
      integer, intent(in) :: n, d      ! number of variables and dimensions
      integer :: j
      u_index_to_u_offset = 0
      do j=1,d
         u_index_to_u_offset = n*u_index_to_u_offset + idx%coor(j)-1
      end do
      u_index_to_u_offset = u_index_to_u_offset + 1
   end function u_index_to_u_offset

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
         do while (i <= l)
            if (a(i) /= v) then
               i = i+1
            else
               exit
            end if
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
         tmp = u_index_to_f_index(tmp)
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

   ! ! Specialized code for local_state_space_iteration_3
   ! ! Considering the folded tensor gᵥᵥ, for each folded offset,
   ! ! fills (i) the corresponding index, (ii) the corresponding 
   ! ! unfolded offset in the corresponding unfolded tensor
   ! ! and (iii) the number of equivalent unfolded indices the folded index
   ! ! associated with the folded offset represents
   ! subroutine index_2(indices, uoff, neq, q)
   !    integer, intent(in) :: q ! size of v
   !    integer, dimension(:), intent(inout) :: uoff, neq ! list of corresponding unfolded offsets and number of equivalent unfolded indices
   !    type(index), dimension(:), intent(inout) :: indices ! list of folded indices
   !    integer :: m, j
   !    m = q*(q+1)/2 ! total number of folded indices : ⎛q+2-1⎞
   !                                                   ! ⎝  2  ⎠
   !    uoff(1) = 1
   !    neq(1) = 1
   !    ! offsets such that j ∈ { 2, ..., q } are associated with
   !    ! indices (1, α), α ∈ { 2, ..., q }
   !    do j=2,q
   !       neq(j) = 2
   !    end do
   ! end subroutine index_2

   ! In order to list folded indices α = (α₁,…,αₘ) with αᵢ ∈ { 1, ..., n },
   ! at least 2 algorithms exist: a recursive one and an iterative one.
   ! The recursive algorithm list_folded_indices(n,m,q) that returns 
   ! the list of all folded indices α = (α₁,…,αₘ) with αᵢ ∈ { 1+q, ..., n+q } works as follows:
   ! if n=0, return an empty list
   ! else if m=0, return the list containing the sole zero-sized index 
   ! otherwise,  
   ! return the concatenation of ([1+q, ℓ] for ℓ ∈ list_folded_indices(n,m-1,q))
   ! and list_folded_indices(n-1,m,1,q+1)]
   ! A call to list_folded_indices(n,m,0) then returns the list
   ! of folded indices α = (α₁,…,αₘ) with αᵢ ∈ { 1, ..., n }
   ! The problem with recursive functions is that the compiler may manage poorly
   ! the stack, which slows down the function's execution
   ! recursive function list_folded_indices(n, m, q) result(list)
   !    integer :: n, m, q 
   !    type(index), allocatable, dimension(:) :: list, temp
   !    integer :: j
   !    if (m==0) then
   !       list = [index(0)]
   !    elseif (n == 0) then
   !       allocate(list(0))
   !    else
   !       temp = list_folded_indices(n,m-1,q)
   !       list = [(index([1+q,temp(j)%coor]), j=1, size(temp)), list_folded_indices(n-1,m,q+1)]
   !    end if
   ! end function list_folded_indices

   ! Considering the folded tensor gᵥᵐ, for each folded offset,
   ! fills the lists of (i) the corresponding index, (ii) the corresponding 
   ! unfolded offset in the corresponding unfolded tensor
   ! and (iii) the number of equivalent unfolded indices the folded index
   ! (associated with the folded offset) represents
   ! The algorithm to get the folded index associated with a folded offset 
   ! relies on the definition of the lexicographic order. 
   ! Considering α = (α₁,…,αₘ) with αᵢ ∈ { 1, ..., n },
   ! the next index α' is such that there exists i that verifies
   ! αⱼ = αⱼ' for all j < i, αᵢ' > αᵢ. Note that all the coordinates
   ! αᵢ', ... , αₘ' need to be as small as the lexicographic order allows
   ! for α' to immediately follow α.
   ! Suppose j is the latest incremented coordinate: 
   ! if αⱼ < n, then αⱼ' = αⱼ + 1
   ! otherwise αⱼ = n, set αₖ' =  αⱼ₋₁ + 1 for all k ≥ j-1
   ! if αⱼ₋₁ = n, set j := j-1  
   ! otherwise, set j := m
   ! The algorithm to count the number of equivalent unfolded indices
   ! works as follows. A folded index can be written as α = (x₁, ..., x₁, ..., xₚ, ..., xₚ)
   ! such that x₁ < x₂ < ... < xₚ. Denote kᵢ the number of coordinates equal to xᵢ.
   ! The number of unfolded indices equivalent to α is c(α) = ⎛         m       ⎞
   !                                                          ⎝ k₁, k₂, ..., kₚ ⎠
   ! Suppose j is the latest incremented coordinate.
   ! If αⱼ < n, then αⱼ' = αⱼ + 1, k(αⱼ) := k(αⱼ)-1, k(αⱼ') := k(αⱼ')+1.
   ! In this case, c(α') = c(α)*(k(αⱼ)+1)/k(αⱼ')
   ! otherwise, αⱼ = n: set αₖ' =  αⱼ₋₁ + 1 for all k ≥ j-1,
   ! k(αⱼ₋₁) := k(αⱼ₋₁)-1, k(n) := 0, k(αⱼ₋₁') = m-(j-1)+1
   ! In this case, we compute c(α') with the multinomial formula above
   ! Finally, the algorithm that returns the unfolded offset of a given folded index works
   ! as follows. Suppose j is the latest incremented coordinate and off(α) is the unfolded offset
   ! associated with index α:
   ! if αⱼ < n, then αⱼ' = αⱼ + 1 and off(α') = off(α)+1
   ! otherwise, αⱼ = n: set αₖ' =  αⱼ₋₁ + 1 for all k ≥ j-1
   ! and off(α') can be computed using the u_index_to_u_offset routine
   subroutine folded_offset_loop(ind, nbeq, off, n, m, p)
      type(index), dimension(:), intent(inout) :: ind ! list of indices
      integer, dimension(:), intent(inout) :: nbeq, off ! lists of numbers of equivalent indices and of offsets
      integer, intent(in) :: n, m
      type(pascal_triangle), intent(in) :: p
      integer :: j, lastinc, k(n)
      ind(1) = index(m, 1)
      nbeq(1) = 1
      k = 0
      k(1) = m
      off(1) = 1
      j = 2
      lastinc = m
      do while (j <= size(ind))
         ind(j) = index(ind(j-1)%coor)
         if (ind(j-1)%coor(lastinc) == n) then
            ind(j)%coor(lastinc-1:m) = ind(j-1)%coor(lastinc-1)+1
            k(ind(j-1)%coor(lastinc-1)) = k(ind(j-1)%coor(lastinc-1))-1
            k(n) = 0
            k(ind(j)%coor(lastinc-1)) = m - (lastinc-1) + 1
            nbeq(j) = multinomial(k,m,p)
            off(j) = u_index_to_u_offset(ind(j), n, m)
            if (ind(j)%coor(m) == n) then 
               lastinc = lastinc-1
            else
               lastinc = m
            end if
         else
            ind(j)%coor(lastinc) = ind(j-1)%coor(lastinc)+1
            k(ind(j)%coor(lastinc)) = k(ind(j)%coor(lastinc))+1
            nbeq(j) = nbeq(j-1)*k(ind(j-1)%coor(lastinc))/k(ind(j)%coor(lastinc))
            k(ind(j-1)%coor(lastinc)) = k(ind(j-1)%coor(lastinc))-1
            off(j) = off(j-1)+1
         end if
         j = j+1
      end do
   end subroutine folded_offset_loop 

   ! The following routine computes the partitions of all integers up to integer
   ! n with number of parts up to n, where n is the size of the input
   ! partition_triangle leading partitions. The partitions are stored in the
   ! lexicographic order. Suppose we want to partition the
   ! integer n using k parts, i.e. solve the problem 𝓟(k,n). Two cases arise:
   ! (i) the partition begins with a one, i.e. is of the form (1, α₁, ...,
   ! αₖ₋₁). In this case, (α₁, ..., αₖ₋₁) is a partition of n-1 with k-1 parts,
   ! i.e. solves the problem 𝓟(k-1,n-1). (ii) Otherwise, if (α₁, ..., αₖ) is a
   ! partition of integer n using k parts, then (α₁-1, ..., αₖ-1) is a partition
   ! of integer n-k using k parts, i.e. solves the problem 𝓟(k,n-k). In other
   ! words, solutions to the problem 𝓟(k,n) are (a) solutions to the problem
   ! 𝓟(k-1,n-1) with a 1 appended up front, and (b) solutions to the problem
   ! 𝓟(k,n-k) with a 1 added to all its indices. Denoting p(k,n) the cardinal
   ! of 𝓟(k,n), it thus verifies p(k,n)=p(k-1,n-1)+p(k,n-k).
   ! A partition of n with k parts can be written as 
   ! α = (α₁, ..., α₁, ..., αₚ, ..., αₚ) such that α₁ < α₂ < ... < αₚ. 
   ! Denote kᵢ the number of coordinates equal to αᵢ. The
   ! number of partitions that are permuted version of α 
   ! is c(α;k) = ⎛       k         ⎞.
   !             ⎝ k₁, k₂, ..., kₚ ⎠
   ! The partitions generated through (b) represent the same number of permuted
   ! partitions. If α solves 𝓟(k,n-k), α'= α.+1 is such that c(α';k)=c(α;k).
   ! As for (a), two cases arise. If α that solves 𝓟(k-1,n-1) is such that
   ! α₁=1, then α'=[1 α] that solves 𝓟(k,n), then
   ! c(α';k) = ⎛       k       ⎞ = ⎛       k-1+1       ⎞ = k*c(α;k-1)/(k₁+1)
   !           ⎝ k₁', ..., kₚ' ⎠   ⎝ k₁+1, k₂, ..., kₚ ⎠
   ! Otherwise, c(α';k)=k*c(α;k-1)
   subroutine fill_partition_triangle(parts)
      type(partition_triangle), dimension(:,:), intent(inout) :: parts 
      integer :: n, k, p_km1_nm1, p_k_nmk, p_k_n, l
      ! Initialization with n=1 
      parts(1,1)%partition = [index(1,1)]
      parts(1,1)%count = [1]
      do n=2,size(parts,1)
         ! print *, 'n:', n
         ! 𝓟(n,n) unique solution is (1,...,1)
         !                            n times
         parts(n,n)%partition = [index(n,1)] 
         parts(n,n)%count = [1]
         ! 𝓟(1,n) unique solution is (n)
         parts(1,n)%partition = [index(1,n)]
         parts(1,n)%count = [1]
         do k=2,n-1
            ! print *, 'k:', k
            p_km1_nm1 = size(parts(k-1,n-1)%partition)
            if (2*k>n) then
               p_k_nmk = 0
            else   
               p_k_nmk = size(parts(k,n-k)%partition)
            end if
            p_k_n = p_km1_nm1+p_k_nmk
            ! print *, 'p_km1_nm1:', p_km1_nm1
            ! print *, 'p_k_nmk:', p_k_nmk
            ! print *, 'p_k_n:', p_k_n
            allocate(parts(k,n)%partition(p_k_n))
            allocate(parts(k,n)%count(p_k_n))
            ! 1 appended up front to 𝓟(k-1,n-1) solutions 
            do l=1,p_km1_nm1
               ! print *, 'l:', l
               allocate(parts(k,n)%partition(l)%coor(k))
               parts(k,n)%partition(l)%coor(1) = 1
               parts(k,n)%partition(l)%coor(2:k) = parts(k-1,n-1)%partition(l)%coor
               if (parts(k-1,n-1)%partition(l)%coor(1) == 1) then
                  parts(k,n)%count(l) = parts(k-1,n-1)%count(l)*k/(get_prefix_length(parts(k-1,n-1)%partition(l), k-1)+1)
               else
                  parts(k,n)%count(l) = parts(k-1,n-1)%count(l)*k
               end if
            end do
            ! 1 added to all components of 𝓟(k,n-k) solutions
            do l=1,p_k_nmk
               ! print *, 'l:', l
               parts(k,n)%partition(l+p_km1_nm1)%coor = parts(k,n-k)%partition(l)%coor+1
               parts(k,n)%count(l+p_km1_nm1) = parts(k,n-k)%count(l)
            end do
         end do
      end do
   end subroutine fill_partition_triangle

end module partitions

! gfortran -o partitions partitions.f08 pascal.f08 sort.f08
! ./partitions
! program test
!    use partitions
!    use pascal
!    implicit none (type, external)
!    type(index) :: uidx, fidx, i1, i2
!    integer, dimension(:), allocatable :: folded
!    integer :: i, uj, n, d, j, nb_folded_idcs, k, l, m
!    type(pascal_triangle) :: p
!    type(index), dimension(:), allocatable :: list_folded_idcs
!    integer, dimension(:), allocatable :: nbeq, off
!    type(partition_triangle), allocatable, dimension(:,:) :: t

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
!    print '(3i2)', (uidx%coor(i), i=1,d) ! should display 0 2 1

!    ! f_index_to_f_offset
!    fidx = u_index_to_f_index(uidx)
!    print '(i2)', f_index_to_f_offset(fidx, n, d, p) ! should display 5

!    ! /=
!    i1 = index((/1,2,3,4,5/))
!    i2 = index((/1,2,3,4,6/))
!    if (i1 /= i2) then
!       print *, "Same!"
!    else
!       print *, "Different!"
!    end if

!    ! fill_folded_indices
!    ! allocate(folded(n**d))
!    ! call fill_folded_indices(folded,n,d,p)
!    ! print *, "Matching offsets unfolded -> folded"
!    ! print '(1000i4)', (i, i=1,n**d)
!    ! print '(1000i4)', (folded(i), i=1,n**d)

!    n = 3
!    d = 3
!    p = pascal_triangle(n+d-1)
!    nb_folded_idcs = get(d,n+d-1,p)
!    ! recursive list_folded_indices
!    ! list_folded_idcs = list_folded_indices(n, d, 0)
!    ! print '(4i2)', ((list_folded_idcs(i)%coor(j), j=1,d), i=1,nb_folded_idcs)

!    ! iterative list_folded_indices
!    allocate(list_folded_idcs(nb_folded_idcs), nbeq(nb_folded_idcs), off(nb_folded_idcs))
!    call folded_offset_loop(list_folded_idcs, nbeq, off, n, d, p)
!    print '(3i2)', ((list_folded_idcs(i)%coor(j), j=1,d), i=1,nb_folded_idcs)
!    print '(i3)', (nbeq(i), i=1,nb_folded_idcs)
!    print '(i4)', (off(i), i=1,nb_folded_idcs)

!    ! Triangle of integer partitions
!    d = 7
!    allocate(t(d,d))
!    call fill_partition_triangle(t)
!    do n=1,d
!       do k=1,n
!          print *, '(k,n):', k, n
!          do l=1, size(t(k,n)%partition)
!             print '(10i2)', (t(k,n)%partition(l)%coor(m), m=1,k)
!             print '(i2)', t(k,n)%count(l)
!          end do
!       end do
!    end do

! end program test