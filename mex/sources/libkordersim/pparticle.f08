! Copyright © 2021-2022 Dynare Team
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

! Routines and data structures for multithreading over particles in local_state_space_iteration_k
module pparticle
   use iso_c_binding
   use simulation
   use matlab_mex

   implicit none

   type tdata
      integer :: nm, nys, endo_nbr, nvar, order, nrestricted, nparticles 
      real(real64), allocatable :: yhat(:,:), e(:,:), ynext(:,:), ys_reordered(:), restrict_var_list(:)
      type(pol), dimension(:), allocatable :: udr
   end type tdata

   type(tdata) :: thread_data

contains

   subroutine thread_eval(arg) bind(c)
      type(c_ptr), intent(in), value :: arg
      integer, pointer :: im
      integer :: i, j, start, end, q, r, ind
      type(horner), dimension(:), allocatable :: h
      real(real64), dimension(:), allocatable :: dyu

      ! Checking that the thread number got passed as argument
      if (.not. c_associated(arg)) then
         call mexErrMsgTxt("No argument was passed to thread_eval")
      end if
      call c_f_pointer(arg, im)

      ! Allocating local arrays
      allocate(h(0:thread_data%order), dyu(thread_data%nvar)) 
      do i=0, thread_data%order
         allocate(h(i)%c(thread_data%endo_nbr, thread_data%nvar**i))
      end do

      ! Specifying bounds for the curent thread
      q = thread_data%nparticles / thread_data%nm
      r = mod(thread_data%nparticles, thread_data%nm)
      start = (im-1)*q+1
      if (im < thread_data%nm) then
         end = start+q-1
      else
         end = thread_data%nparticles
      end if

      ! Using the Horner algorithm to evaluate the decision rule at the chosen yhat and epsilon
      do j=start,end
         dyu(1:thread_data%nys) = thread_data%yhat(:,j) 
         dyu(thread_data%nys+1:) = thread_data%e(:,j) 
         call eval(h, dyu, thread_data%udr, thread_data%endo_nbr, thread_data%nvar, thread_data%order)
         do i=1,thread_data%nrestricted
            ind = int(thread_data%restrict_var_list(i))
            thread_data%ynext(i,j) = h(0)%c(ind,1) + thread_data%ys_reordered(ind)
         end do
      end do

   end subroutine thread_eval

end module pparticle
