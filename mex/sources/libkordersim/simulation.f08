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
   use partitions
   use lapack
   implicit none (type, external)

   ! Used to store the folded decision rule tensors
   type :: pol
      real(real64), allocatable :: g(:,:)
   end type pol

   ! Used to store the different contracted tensors used in the Horner algorithm
   ! type :: horner
   !    real(real64), pointer, contiguous :: c(:,:), c_flat(:)
   ! end type horner
   type :: horner
      real(real64), allocatable :: c(:,:)
   end type horner


contains

   ! ! With MATMUL
   ! ! Horner evaluation of the polynomial with coefficients stored in udr at the point dyu
   ! subroutine eval(h, dyu, udr, ny, nvar, order)
   !    type(horner), dimension(0:order), intent(inout) :: h
   !    real(real64), dimension(nvar), intent(in) :: dyu
   !    type(pol), intent(in) :: udr(0:order)
   !    integer, intent(in) :: ny, nvar, order
   !    integer :: d
   !    if (order == 1) then
   !       h(1)%c = udr(1)%g
   !    else
   !       call contract(h(order-1)%c, udr(order)%g, dyu, ny, nvar, order)
   !    end if
   !    do d=order-1,1,-1
   !       if (d > 1) then
   !          h(d)%c = h(d)%c + udr(d)%g
   !          call contract(h(d-1)%c, h(d)%c, dyu, ny, nvar, d)
   !       else
   !          h(d)%c = h(d)%c + udr(1)%g
   !       end if
   !    end do
   !    h(0)%c(:,1) = matmul(h(1)%c, dyu) + udr(0)%g(:,1)
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
      type(horner), dimension(0:order), intent(inout) :: h
      real(real64), dimension(nvar), intent(in) :: dyu
      type(pol), intent(in) :: udr(0:order)
      integer, intent(in) :: ny, nvar, order
      integer :: d
      if (order == 1) then
         h(1)%c = udr(1)%g
      else
         call contract(h(order-1)%c, udr(order)%g, dyu, ny, nvar, order)
      end if
      do d=order-1,1,-1
         if (d > 1) then
            h(d)%c = h(d)%c + udr(d)%g
            call contract(h(d-1)%c, h(d)%c, dyu, ny, nvar, d)
         else
            h(d)%c = h(d)%c + udr(1)%g
         end if
      end do
      h(0)%c = udr(0)%g
      call mult_vec(1.0_real64, h(1)%c, dyu, 1.0_real64, h(0)%c(:,1))
   end subroutine eval

   ! Contracts the unfolded tensor t with respect to the vector x and stores the result in c
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

end module simulation

