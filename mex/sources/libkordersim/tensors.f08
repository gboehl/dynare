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

module tensors
   use iso_c_binding
   use partitions
   ! use matlab_mat
   implicit none (type, external)

   ! A type to contain a folded or unfolded tensor [gᵥᵐ]
   type tensor
      real(real64), pointer, contiguous, dimension(:,:) :: m
      ! real(real64), dimension(:,:), allocatable :: m
   end type tensor 

   ! A type to contain the unfolded tensor [gᵥᵐ]. For each offset j, ind(j)
   ! contains the associated unfolded index and s(j) stores the number of state 
   ! variables in the index, i.e. if v = (x,u[,σ]) where x is the s-sized state 
   ! variables vector and ind(j) = (α₁,…,αₘ), s(j) stores the number of 
   ! coordinates αᵢ such that  αᵢ ≤ s.
   type unfolded_tensor
      type(index), dimension(:), allocatable :: ind
      integer, dimension(:), allocatable :: s
      real(real64), dimension(:,:), pointer, contiguous :: m
   end type unfolded_tensor

   ! A type to contain the folded tensor [gᵥᵐ]. For each offset j, ind(j)
   ! contains the associated folded index, count(j) contains the number of
   ! equivalent unfolded indices and s(j) stores the number of state variables
   ! in the index, i.e. if v = (x,u[,σ]) where x is the s-sized state variables
   ! vector and ind(j) = (α₁,…,αₘ), s(j) stores the number of coordinates αᵢ
   ! such that  αᵢ ≤ s.
   type folded_tensor
      type(index), dimension(:), allocatable :: ind
      integer, dimension(:), allocatable :: count
      integer, dimension(:), allocatable :: s
      real(real64), dimension(:,:), pointer, contiguous :: m
   end type folded_tensor

contains

   ! Fills the ind and s members of a list of unfolded tensors g
   ! with a number of state variables ns. More precisely, g contains the
   ! information associated with the unfolded tensors [gᵥᵐ] with m=1,...,q
   ! 1. Filling g%ind
   ! The unfolded indices of [gᵥᵐ⁺¹] are the unfolded indices of [gᵥᵐ]
   ! with all numbers 1,...,n appended upfront. In other words, any
   ! index of [gᵥᵐ⁺¹] writes (α₁,α) with α₁∈ { 1, ..., n } and α 
   ! an index of [gᵥᵐ]
   ! 2. Filling g%s
   ! s(α₁,α) is s(α)+1 if α₁ ≤ ns and s(α) otherwise 
   subroutine fill_list_unfolded_tensor(g, s, n)
      type(unfolded_tensor), dimension(:), intent(inout) :: g
      integer, intent(in) :: s, n ! number of state variables and size of v
      integer :: q, m, j, k, l, length
      q = size(g)
      ! Initialization
      m = 1
      length = n
      allocate(g(1)%ind(length), g(1)%s(length))
      do j=1,n
         g(1)%ind(j) = index(1,j)
         if (j<=s) then
            g(1)%s(j) = 1
         else
            g(1)%s(j) = 0
         end if
      end do
      ! Transmission
      do m=2,q
         length = n*length
         allocate(g(m)%ind(length), g(m)%s(length))
         l = 1
         do j=1,n
            do k=1,size(g(m-1)%ind)
               g(m)%ind(l) = index(m)
               g(m)%ind(l)%coor(1) = j
               g(m)%ind(l)%coor(2:m) = g(m-1)%ind(k)%coor
               if (j<=s) then
                  g(m)%s(l) = g(m-1)%s(k)+1
               else
                  g(m)%s(l) = g(m-1)%s(k)
               end if
               l = l+1
            end do
         end do
      end do
   end subroutine fill_list_unfolded_tensor

   ! Allocates the members ind, count and s of a tensor [gᵥᵐ] denoted g 
   ! using the Pascal triangle p. n is the size of v.
   subroutine allocate_folded_tensor(g, m, n, p)
      type(folded_tensor), intent(inout) :: g
      integer, intent(in) :: m, n ! 
      type(pascal_triangle),intent(in) :: p
      integer :: d
      d = get(m,n+m-1,p)
      allocate(g%ind(d), g%count(d), g%s(d))
   end subroutine allocate_folded_tensor
   
   ! Fills the already allocated members ind, count and s of a tensor [gᵥᵐ] 
   ! denoted g using the Pascal triangle p. n is the size of v.
   ! 1. Filling g%ind
   ! The algorithm to get the folded index associated with a folded offset 
   ! relies on the definition of the lexicographic order. 
   ! Consideri appended upfront.ng α = (α₁,…,αₘ) with αᵢ ∈ { 1, ..., n },
   ! the next index α' is such that there exists i that verifies
   ! αⱼ = αⱼ' for all j < i, αᵢ' > αᵢ. Note that all the coordinates
   ! αᵢ', ... , αₘ' need to be as small as the lexicographic order allows
   ! for α' to immediately follow α.
   ! Suppose j is the latest incremented coordinate: 
   ! if αⱼ < n, then αⱼ' = αⱼ + 1
   ! otherwise αⱼ = n, set αₖ' =  αⱼ₋₁ + 1 for all k ≥ j-1
   ! if αⱼ₋₁ = n, set j := j-1  
   ! otherwise, set j := m
   ! 2. Filling g%count
   ! The algorithm to count the number of equivalent unfolded indices
   ! works as follows. A folded index can be written as 
   ! α = (x₁, ..., x₁, ..., xₚ, ..., xₚ) such that x₁ < x₂ < ... < xₚ.
   ! Denote kᵢ the number of coordinates equal to xᵢ.
   ! The number of unfolded indices equivalent to α is 
   ! c(α) = ⎛        m        ⎞
   !        ⎝ k₁, k₂, ..., kₚ ⎠
   ! k is an array such that k(ℓ,j) contains the number of coordinates
   ! equal to ℓ for the folded index asociated with offset j. 
   ! Suppose j is the latest incremented coordinate.
   ! If αⱼ < n, then αⱼ' = αⱼ + 1, k(αⱼ) := k(αⱼ)-1, k(αⱼ') := k(αⱼ')+1.
   ! In this case, c(α') = c(α)*(k(αⱼ)+1)/k(αⱼ')
   ! otherwise, αⱼ = n: set αₖ' =  αⱼ₋₁ + 1 for all k ≥ j-1,
   ! k(αⱼ₋₁) := k(αⱼ₋₁)-1, k(n) := 0, k(αⱼ₋₁') = m-(j-1)+1
   ! In this case, we compute c(α') with the multinomial formula above
   ! 3. Filling g%s
   ! Suppose j is the latest incremented coordinate.
   ! If αⱼ < n, then αⱼ' = αⱼ + 1: s' = s-1 if αⱼ = ns and s'=s otherwise
   ! Otherwise, αⱼ = n: set αₖ' =  αⱼ₋₁ + 1 for all k ≥ j-1. Thus, 
   ! s' = m if αⱼ₋₁ < ns; s'=s-1 if αⱼ₋₁ = ns; s'=s  otherwise
   subroutine fill_folded_tensor(g, k, m, s, n, p)
      type(folded_tensor), intent(inout) :: g
      integer, contiguous, intent(inout) :: k(:,:)
      integer, intent(in) :: m, s, n
      type(pascal_triangle) :: p
      integer :: j, lastinc
      g%ind(1) = index(m, 1)
      g%count(1) = 1
      g%s(1) = m
      k = 0
      k(:,1) = 0
      k(1,1) = m
      j = 2
      lastinc = m
      do while (j <= size(g%ind))
         g%ind(j) = index(g%ind(j-1)%coor)
         k(:,j) = k(:,j-1)
         if (g%ind(j-1)%coor(lastinc) == n) then
            g%ind(j)%coor(lastinc-1:m) = g%ind(j-1)%coor(lastinc-1)+1
            k(g%ind(j-1)%coor(lastinc-1),j) = k(g%ind(j-1)%coor(lastinc-1),j-1)-1
            k(n,j) = 0
            k(g%ind(j)%coor(lastinc-1),j) = m - (lastinc-1) + 1
            g%count(j) = multinomial(k(:,j),m,p)
            if (g%ind(j-1)%coor(lastinc-1) < s) then
               g%s(j) = m 
            elseif (g%ind(j-1)%coor(lastinc-1) == s) then
               g%s(j) = g%s(j-1)-1
            else
               g%s(j) = g%s(j-1)
            end if
            if (g%ind(j)%coor(m) == n) then 
               lastinc = lastinc-1
            else
               lastinc = m
            end if
         else
            g%ind(j)%coor(lastinc) = g%ind(j-1)%coor(lastinc)+1
            k(g%ind(j)%coor(lastinc),j) = k(g%ind(j)%coor(lastinc),j-1)+1
            g%count(j) = g%count(j-1)*k(g%ind(j-1)%coor(lastinc),j)/k(g%ind(j)%coor(lastinc),j)
            k(g%ind(j-1)%coor(lastinc),j) = k(g%ind(j-1)%coor(lastinc),j-1)-1
            if (g%ind(j-1)%coor(lastinc) == s) then
               g%s(j) = g%s(j-1)-1 
            else
               g%s(j) = g%s(j-1)
            end if
         end if
         j = j+1
      end do
   end subroutine

   ! Fills a tensor [gᵥᵐ] given the tensor [gᵥᵐ⁺¹] 
   ! 1. Filling g%ind 
   ! If (α₁,…,αₘ₊₁) is a folded index for [gᵥᵐ⁺¹], then (α₂,…,αₘ₊₁) is a folded
   ! index of [gᵥᵐ]. The folded indices of [gᵥᵐ] are thus the tails of the first
   ! ⎛ m+n-1 ⎞ folded indices of [gᵥᵐ⁺¹]
   ! ⎝   m   ⎠
   ! 2. Filling g%count
   ! If α=(α₁,…,αₘ₊₁) is a folded index for [gᵥᵐ⁺¹], then α'=(α₂,…,αₘ₊₁) is a
   ! folded index of [gᵥᵐ]. We thus have c(α') = c(α)*k(α₁)/(m+1) and 
   ! perform k(α₁):=k(α₁)-1.
   ! 3. Filling g%s
   ! If α=(α₁,…,αₘ₊₁) is a folded index for [gᵥᵐ⁺¹], then α'=(α₂,…,αₘ₊₁) is a
   ! folded index of [gᵥᵐ]. We thus have s(α') = s(α)-1 if α₁ ≤ s and 
   ! s(α') = s(α) otherwise
   subroutine fill_folded_tensor_backward(g_m, k, g_mp1, m, s)
      type(folded_tensor), intent(inout) :: g_m ! tensor [gᵥᵐ]
      ! k array from previous fill_folded_tensor_backward or 
      ! fill_folded_tensor call. See definition in fill_folded_tensor
      integer, contiguous, intent(inout) :: k(:,:)
      type(folded_tensor), intent(in) :: g_mp1 ! tensor [gᵥᵐ⁺¹]
      integer, intent(in) :: m, s ! s is the number of state variables
      integer :: j
      do j=1,size(g_m%ind)
         g_m%ind(j)%coor = g_mp1%ind(j)%coor(2:m+1)
         g_m%count(j) = g_mp1%count(j)*k(g_mp1%ind(j)%coor(1),j)/(m+1)
         k(g_mp1%ind(j)%coor(1),j) = k(g_mp1%ind(j)%coor(1),j)-1
         if (g_mp1%ind(j)%coor(1) <= s) then
            g_m%s(j) = g_mp1%s(j)-1
         else
            g_m%s(j) = g_mp1%s(j)
         end if
      end do
   end subroutine

   ! Fills the ind, count and s members of a list of folded tensors g
   ! with a number of state variables ns. More precisely, g contains the
   ! information associated with the folded tensors [gᵥᵐ] with m=1,...,q
   subroutine fill_list_folded_tensor(g, s, n, p)
      type(folded_tensor), dimension(:), intent(inout) :: g
      integer, intent(in) :: s, n ! number of state variables and size of v
      type(pascal_triangle) :: p
      integer :: q, m 
      integer, dimension(:,:), allocatable :: k
      q = size(g)
      ! Case m = q
      m = q
      call allocate_folded_tensor(g(m), m, n, p)
      allocate(k(n,size(g(m)%ind)))
      call fill_folded_tensor(g(m), k, m, s, n, p) 
      ! Case m < q
      do m=q-1,1,-1
         call allocate_folded_tensor(g(m), m, n, p)
         call fill_folded_tensor_backward(g(m), k, g(m+1), m, s)
      end do
   end subroutine fill_list_folded_tensor
end module tensors

! After putting a comment on MATLAB-related lines
! gfortran -o tensors tensors.f08 sort.f08 pascal.f08 partitions.f08
! ./tensors
! program test
!    use pascal
!    use tensors
!    implicit none (type, external)
!    type(pascal_triangle) :: p
!    type(folded_tensor) :: g, h
!    integer :: n, m, s, i, j
!    integer, allocatable, dimension(:,:) :: k
 
!    n = 3
!    m = 3
!    s = 2
!    p = pascal_triangle(n+m-1)

!    call allocate_folded_tensor(g, m, n, p)
!    allocate(k(n, size(g%ind)))
!    call fill_folded_tensor(g, k, c_null_ptr, m, s, n, p)

!    ! List of folded indices, counts of equivalent unfolded indices
!    ! and counts of state variables of [gᵥᵐ]
!    ! Folded indices and offsets
!    ! 0,0,0  1     1,1,1   7      2,2,2   10
!    ! 0,0,1  2     1,1,2   8
!    ! 0,0,2  3     1,2,2   9
!    ! 0,1,1  4
!    ! 0,1,2  5
!    ! 0,2,2  6
!    do i=1, size(k,2)
!       print '(3i2)', (g%ind(i)%coor(j), j=1,m)
!       print '(i2)', g%count(i)
!       print '(i2)', g%s(i)
!    end do

!    call allocate_folded_tensor(h, m-1, n, p)
!    call fill_folded_tensor_backward(h, k, g, c_null_ptr, m-1, s)
!    ! List of folded indices, counts of equivalent unfolded indices
!    ! and counts of state variables of [gᵥᵐ⁻¹]
!    do i=1, size(h%ind)
!       print '(2i2)', (h%ind(i)%coor(j), j=1,m-1)
!       print '(i2)', h%count(i)
!       print '(i2)', h%s(i)
!    end do

! end program test