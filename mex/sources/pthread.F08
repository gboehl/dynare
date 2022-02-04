! Copyright © 2022 Dynare Team
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

! Wrapper around pthread_create and pthread_join C routines for POSIX multithreading
module pthread
   use iso_c_binding
   implicit none
   private

   public :: c_pthread_create
   public :: c_pthread_join

   ! SIZEOF_PTHREAD_T is set by the AC_CHECK_SIZEOF routine
   integer, parameter :: pthread_size = SIZEOF_PTHREAD_T

   type, bind(c), public :: c_pthread_t
      private
      character(kind=c_char) :: hidden(pthread_size)
   end type c_pthread_t

   interface

      function c_pthread_create(thread, attr, start_routine, arg) bind(c, name='pthread_create')
         import :: c_int, c_ptr, c_funptr, c_pthread_t
         implicit none
         type(c_pthread_t), intent(inout)     :: thread
         type(c_ptr), intent(in), value :: attr
         type(c_funptr), intent(in), value :: start_routine
         type(c_ptr), intent(in), value :: arg
         integer(kind=c_int) :: c_pthread_create
      end function c_pthread_create

      function c_pthread_join(thread, value_ptr) bind(c, name='pthread_join')
         import :: c_int, c_ptr, c_pthread_t
         implicit none
         type(c_pthread_t), intent(in), value :: thread
         type(c_ptr),       intent(in) :: value_ptr
         integer(kind=c_int) :: c_pthread_join
      end function c_pthread_join

    end interface
end module pthread
