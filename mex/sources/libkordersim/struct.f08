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
module struct
   use iso_fortran_env
   use iso_c_binding
   use matlab_mex
   implicit none

   contains

      type(integer) function get_int_field(struct, field)
         type(c_ptr), intent(in) :: struct
         character(*), intent(in) :: field
         type(c_ptr) :: tmp
         tmp = mxGetField(struct, 1_mwIndex, field)
         if (.not. (c_associated(tmp) .and. mxIsScalar(tmp) .and. mxIsNumeric(tmp))) then
            call mexErrMsgTxt("Field "//field//" should be a numeric scalar")
         end if
         get_int_field = int(mxGetScalar(tmp))
      end function get_int_field

end module struct