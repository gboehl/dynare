! Necessary routines and functions to carry out simulations
!
! A first step is to get the associated 

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

module simulation
   use iso_fortran_env
   use tensors
   use blas
   use matlab_mex
   implicit none (type, external)

   ! Used to store the simulated pruned state space
   type pruned
      real(real64), allocatable :: y(:,:)
   end type pruned

contains

   ! ! With MATMUL
   ! ! Horner evaluation of the polynomial with coefficients stored in udr at the point dyu
   ! subroutine eval(h, dyu, udr, ny, nvar, order)
   !    type(tensor), dimension(0:order), intent(inout) :: h
   !    real(real64), dimension(nvar), intent(in) :: dyu
   !    type(tensor), intent(in) :: udr(0:order)
   !    integer, intent(in) :: ny, nvar, order
   !    integer :: d
   !    if (order == 1) then
   !       h(1)%c = udr(1)%m
   !    else
   !       call contract(h(order-1)%c, udr(order)%g, dyu, ny, nvar, order)
   !    end if
   !    do d=order-1,1,-1
   !       if (d > 1) then
   !          h(d)%c = h(d)%c + udr(d)%m
   !          call contract(h(d-1)%c, h(d)%c, dyu, ny, nvar, d)
   !       else
   !          h(d)%c = h(d)%c + udr(1)%m
   !       end if
   !    end do
   !    h(0)%c(:,1) = matmul(h(1)%c, dyu) + udr(0)%m(:,1)
   ! end subroutine eval

   ! ! Contracts the unfolded tensor t with respect to the vector x and stores the result in c
   ! subroutine contract(c, t, x, nrows, nvar, d)
   !    integer, intent(in) :: nrows, nvar, d
   !    real(real64), dimension(nrows, nvar**d), intent(in) :: t
   !    real(real64), dimension(nvar), intent(in) :: x
   !    real(real64), dimension(nrows, nvar**(d-1)), intent(inout) :: c(:,:)
   !    real(real64), dimension(nrows) :: tmp
   !    integer :: i
   !    do i=1,nvar**(d-1)
   !       tmp = matmul(t(:,(i-1)*nvar+1:i*nvar), x)
   !       c(:,i) = tmp
   !    end do
   ! end subroutine contract

   ! With DGEMV   
   ! Modifies y such that y := alpha*A*x + beta*y, where alpha, beta are scalars,
   ! x and y are vectors, A is a m-by-n matrix
   subroutine mult_vec(alpha, A, x, beta, y)
      real(real64), intent(in) :: alpha, beta, x(:), A(:,:)
      real(real64), intent(inout) :: y(:)
      call dgemv("N", int(size(A,1), blint), int(size(A,2), blint), alpha, A, int(size(A,1), blint), x, 1_blint, beta, y, 1_blint)
   end subroutine mult_vec

   ! Horner evaluation of the polynomial with coefficients stored in udr at the point dyu
   subroutine eval(h, dyu, udr, ny, nvar, order)
      type(tensor), dimension(0:order), intent(inout) :: h
      real(real64), dimension(nvar), intent(in) :: dyu
      type(tensor), intent(in) :: udr(0:order)
      integer, intent(in) :: ny, nvar, order
      integer :: d
      if (order == 1) then
         h(1)%m = udr(1)%m
      else
         call contract(h(order-1)%m, udr(order)%m, dyu, ny, nvar, order)
      end if
      do d=order-1,1,-1
         if (d > 1) then
            h(d)%m = h(d)%m + udr(d)%m
            call contract(h(d-1)%m, h(d)%m, dyu, ny, nvar, d)
         else
            h(d)%m = h(d)%m + udr(1)%m
         end if
      end do
      h(0)%m = udr(0)%m
      call mult_vec(1.0_real64, h(1)%m, dyu, 1.0_real64, h(0)%m(:,1))
   end subroutine eval

   ! Contracts the unfolded tensor t with respect to the vector x and stores the
   ! result in c.
   subroutine contract(c, t, x, nrows, nvar, d)
      integer, intent(in) :: nrows, nvar, d
      real(real64), dimension(nrows, nvar**d), intent(in) :: t
      real(real64), dimension(nvar), intent(in) :: x
      real(real64), dimension(nrows, nvar**(d-1)), intent(inout) :: c
      integer :: i
      do i=1,nvar**(d-1)
         call mult_vec(1.0_real64, t(:,(i-1)*nvar+1:i*nvar), x, 0.0_real64, c(:,i)) 
      end do
   end subroutine contract

   ! Simulates the model around the steady state using the pruned state space
   subroutine simulate_pruning(sim, dr_mx, ysteady, dy, shocks, order, nstatic, nvar)
      real(real64), dimension(:,:), contiguous, intent(inout) :: sim
      type(c_ptr), intent(in) :: dr_mx
      real(real64), contiguous, intent(in) :: ysteady(:), dy(:) 
      real(real64), contiguous, target, intent(in) :: shocks(:,:)
      integer, intent(in) :: order, nstatic, nvar

      real(real64), dimension(:,:), allocatable :: phat
      real(real64), dimension(:), contiguous, pointer :: u
      type(tensor), dimension(:), allocatable :: psim
      type(partition_triangle), dimension(:,:), allocatable, target :: part
      real(real64) :: prod, sum_part
      type(partition_triangle), pointer :: part_set
      character(kind=c_char, len=10) :: fieldname
      type(c_ptr) :: g_q
      type(unfolded_tensor), dimension(:), allocatable :: g
      integer :: q, l, j, m, r, t, c, endo_nbr, nys, nper
      ! Import the MATLAB decision rule matrices
      q = 1
      allocate(g(1:order))
      do while (q <= order)
         write (fieldname, '(a2, i1)') "g_", q
         g_q = mxGetField(dr_mx, 1_mwIndex, trim(fieldname))
         if (.not. (c_associated(g_q) .and. mxIsDouble(g_q) .and. .not. mxIsComplex(g_q) .and. .not. mxIsSparse(g_q))) then
            call mexErrMsgTxt(trim(fieldname)//" is not allocated in dr.pruning")
         end if
         g(q)%m(1:mxGetM(g_q),1:mxGetN(g_q)) => mxGetPr(g_q)
         q = q+1
      end do
      ! Initialize useful variables
      nys = size(dy)
      endo_nbr = size(ysteady)
      nper = size(sim,2)
      allocate(psim(order), phat(nys,order), part(order,order))
      allocate(psim(1)%m(endo_nbr, 1))
      phat(:,1) = dy
      do l=2,order
         allocate(psim(l)%m(endo_nbr, 1))
         psim(l)%m(:,1) = 0._real64
         phat(:,l) = 0._real64
      end do
      call fill_list_unfolded_tensor(g, nys, nvar+1)
      call fill_partition_triangle(part)
      ! Loop over periods
      do t=2,nper
         u => shocks(:,t-1)
         do l=1,order
            psim(l)%m(:,1) = 0._real64
         end do
         ! Go over the [gᵥˡ] folded tensors
         do l=1,order
            ! Go over the offsets of the matrix spanning [gᵥˡ]
            do j=1,size(g(l)%ind)
               if (g(l)%s(j) > 0) then
                  ! Compute the contribution of g_l-column-multiplied
                  ! terms to x_q when the index of the g_l columns
                  ! involves at least one state variable
                  ! Contributions are non-zero if q ≥ l
                  do q=l,order
                     ! 1. Get the sum-over-integer-partition terms
                     sum_part = 0._real64
                     part_set => part(g(l)%s(j),q-(l-g(l)%s(j)))
                     do r=1,size(part_set%partition)
                        prod = real(part_set%count(r),real64)
                        c = 1
                        do m=1,l
                           if (g(l)%ind(j)%coor(m) <= nys) then
                              prod = prod*phat(g(l)%ind(j)%coor(m),&
                              &part_set%partition(r)%coor(c))
                              c = c+1
                           end if
                           if ((g(l)%ind(j)%coor(m) > nys) .and. (g(l)%ind(j)%coor(m) <= nvar)) then
                              prod = prod*u(g(l)%ind(j)%coor(m)-nys)
                           end if
                        end do
                        sum_part = sum_part+prod
                     end do
                     ! 2. Multiply the sum-over-integer-partition terms
                     ! by the considered column [gᵥˡ] 
                     psim(q)%m(:,1) = psim(q)%m(:,1)+g(l)%m(:,j)*sum_part
                  end do
               else
                  ! Compute the contribution of g_l-column-multiplied
                  ! terms to x_q when the index of the g_l columns
                  ! involves no state variables. Such a contribution only 
                  ! exists when q=l.
                  prod = 1._real64
                  do m=1,l
                     if ((g(l)%ind(j)%coor(m) > nys) .and. (g(l)%ind(j)%coor(m) <= nvar)) then
                        prod = prod*u(g(l)%ind(j)%coor(m)-nys)
                     end if
                  end do
                  psim(l)%m(:,1) = psim(l)%m(:,1)+g(l)%m(:,j)*prod
               end if
            end do
         end do
         ! Add l-order terms to the output simulation for period t
         ! for 1 ≤ l ≤ order
         sim(:,t) = ysteady
         ! Prepare next iteration
         do l=1,order
            phat(:,l) = psim(l)%m(nstatic+1:nstatic+nys,1)
            sim(:,t) = sim(:,t)+psim(l)%m(:,1)
         end do
      end do
   end subroutine simulate_pruning

   ! Simulates the a model around the steady state
   subroutine simulate(sim, dr_mx, ysteady, dy, shocks, order, nstatic, nvar)
      real(real64), dimension(:,:), contiguous, intent(inout) :: sim
      type(c_ptr), intent(in) :: dr_mx
      real(real64), contiguous, intent(in) :: ysteady(:), dy(:), shocks(:,:)
      integer, intent(in) :: order, nstatic, nvar
      character(kind=c_char, len=10) :: fieldname
      type(c_ptr) :: g_d
      type(tensor), dimension(:), allocatable :: fg, ug, h
      integer :: d, endo_nbr, nys, t
      type(uf_matching), dimension(:), allocatable :: matching
      real(real64), dimension(:), allocatable :: dyu
      type(pascal_triangle) :: p
      ! Import the MATLAB decision rule matrices
      d = 0
      allocate(fg(0:order))
      do while (d <= order)
         write (fieldname, '(a2, i1)') "g_", d
         g_d = mxGetField(dr_mx, 1_mwIndex, trim(fieldname))
         if (.not. (c_associated(g_d) .and. mxIsDouble(g_d) .and. .not. mxIsComplex(g_d) .and. .not. mxIsSparse(g_d))) then
            call mexErrMsgTxt(trim(fieldname)//" is not allocated in dr")
         end if
         fg(d)%m(1:mxGetM(g_d),1:mxGetN(g_d)) => mxGetPr(g_d)
         d = d+1
      end do
      ! Put decision rules matrices in the unfolded form
      ! Pinpointing the corresponding offsets between folded and unfolded
      ! tensors
      allocate(h(0:order), ug(0:order), matching(2:order))
      endo_nbr = size(ysteady)
      do d=0,1
         allocate(h(d)%m(endo_nbr,nvar**d))
         ug(d)%m => fg(d)%m
      end do
      ! Compute the useful binomial coefficients from Pascal's triangle
      p = pascal_triangle(nvar+order-1)
      do d=2,order
         allocate(h(d)%m(endo_nbr,nvar**d), ug(d)%m(endo_nbr, nvar**d), matching(d)%folded(nvar**d))
         call fill_folded_indices(matching(d)%folded, nvar, d, p)
         ug(d)%m = fg(d)%m(:,matching(d)%folded)
      end do
      allocate(dyu(nvar))
      ! Getting the predetermined part of the endogenous variable vector 
      nys = size(dy)
      dyu(1:nys) = dy 
      ! Carrying out the simulation
      do t=2,size(sim,2)
         dyu(nys+1:) = shocks(:,t-1)
         ! Using the Horner algorithm to evaluate the decision rule at the 
         ! chosen dyu
         call eval(h, dyu, ug, endo_nbr, nvar, order)
         sim(:,t) = h(0)%m(:,1) + ysteady
         dyu(1:nys) = h(0)%m(nstatic+1:nstatic+nys,1)
      end do
   end subroutine simulate



end module simulation

