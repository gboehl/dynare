! Fortran 2008 interface for a subset of MEX functions
!
! For some functions, exposing the C interface directly is not convenient (e.g.
! when they use C strings, since those need to be null-terminated, or when they
! return pointers to arrays of doubles/chars/integers…). A more Fortran-ish
! interface is provided for those, using a small glue code.
!
! Things to be aware of when adding new interfaces this file:
!
! — The tricky part is to deal with API versioning.
!   • API_VER is for functions which were not versioned before 9.4 (R2018a)
!   • API_VER2 is for functions which were versioned _730 before 9.4 (R2018a)
!     when not using the MX_COMPAT_32 mode
!   For each function, the information can be retrieved from either matrix.h or
!   mex.h from R2014a
! — C passes arguments by value, so the “value” keyword is often needed
! — Strings passed to C must be null terminated (hence a wrapper is needed to
!   append c_null_char)
! — We follow the Fortran convention that indices start at one. Hence, when
!   interfacing a function with a mwIndex argument, it is necessary to write
!   a glue code that substracts one to the index
! — When writing glue code, using the pure C interface as a starting point:
!   • remove the “use” declarations
!   • remove the “value” keywords
!   • convert input character arrays to character(kind=c_char, len=*)
!   • Fortran array pointers returned by the glue code must be marked
!     “contiguous” (which is always the case for C arrays). This will help the
!     Fortran compiler better optimize the code (in some cases, this will avoid
!     array copies)

! Copyright © 2019-2022 Dynare Team
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

#ifdef MATLAB_MEX_FILE
# if MATLAB_VERSION >= 0x0904
#  define API_VER "_800"
#  define API_VER2 "_800"
# else
#  define API_VER ""
#  define API_VER2 "_730"
# endif
#else
! Octave
# define API_VER ""
# define API_VER2 ""
#endif

#include "defines.F08"

!!! C Matrix API
!!! Listed in same order as https://fr.mathworks.com/help/matlab/cc-mx-matrix-library.html
module matlab_mat
  use iso_fortran_env
  use iso_c_binding
  implicit none

  !! C Data Types
  integer, parameter :: mwSize = c_size_t
  integer, parameter :: mwIndex = c_size_t
  integer, parameter :: mwSignedIndex = c_intptr_t
  integer, parameter :: mxLogical = c_bool
  integer, parameter :: mxComplexity = c_int

  integer(mxComplexity), parameter :: mxREAL = 0
  integer(mxComplexity), parameter :: mxCOMPLEX = 1

  interface
     !! mxArray attributes
     logical(c_bool) function mxIsNumeric(pm) bind(c, name="mxIsNumeric"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsNumeric

     logical(c_bool) function mxIsComplex(pm) bind(c, name="mxIsComplex"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsComplex

     integer(c_size_t) function mxGetNumberOfElements(pm) bind(c, name="mxGetNumberOfElements"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetNumberOfElements

     integer(c_size_t) function mxGetM(pm) bind(c, name="mxGetM"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetM

     integer(c_size_t) function mxGetN(pm) bind(c, name="mxGetN"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetN

     logical(c_bool) function mxIsEmpty(pm) bind(c, name="mxIsEmpty"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsEmpty

     !! Create, Query, and Access Data Types

     ! Numeric types
     type(c_ptr) function mxCreateDoubleMatrix(m, n, ComplexFlag) bind(c, name="mxCreateDoubleMatrix"//API_VER2)
       use iso_c_binding
       import :: mwSize, mxComplexity
       integer(mwSize), intent(in), value :: m, n
       integer(mxComplexity), intent(in), value :: ComplexFlag
     end function mxCreateDoubleMatrix

     type(c_ptr) function mxCreateDoubleScalar(value) bind(c, name="mxCreateDoubleScalar"//API_VER)
       use iso_c_binding
       real(c_double), intent(in), value :: value
     end function mxCreateDoubleScalar

     ! Noncomplex Float
     logical(c_bool) function mxIsScalar(array_ptr) bind(c, name="mxIsScalar"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: array_ptr
     end function mxIsScalar

     real(c_double) function mxGetScalar(pm) bind(c, name="mxGetScalar"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetScalar

     logical(c_bool) function mxIsDouble(pm) bind(c, name="mxIsDouble"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsDouble

#if MX_HAS_INTERLEAVED_COMPLEX
     type(c_ptr) function mxGetDoubles_internal(pm) bind(c, name="mxGetDoubles"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetDoubles_internal
#endif

     type(c_ptr) function mxGetPr_internal(pm) bind(c, name="mxGetPr"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetPr_internal

     ! Noncomplex integer
#if MX_HAS_INTERLEAVED_COMPLEX
     type(c_ptr) function mxGetInt32s_internal(pa) bind(c, name="mxGetInt32s"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pa
     end function mxGetInt32s_internal
#endif

     ! Complex Float
#if MX_HAS_INTERLEAVED_COMPLEX
     type(c_ptr) function mxGetComplexDoubles_internal(pa) bind(c, name="mxGetComplexDoubles"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pa
     end function mxGetComplexDoubles_internal
#else
     type(c_ptr) function mxGetPi_internal(pm) bind(c, name="mxGetPi"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetPi_internal
#endif

     ! Sparse
     type(c_ptr) function mxCreateSparse(m, n, nzmax, ComplexFlag) bind(c, name="mxCreateSparse"//API_VER2)
       use iso_c_binding
       import :: mwSize, mxComplexity
       integer(mwSize), intent(in), value :: m, n, nzmax
       integer(mxComplexity), intent(in), value :: ComplexFlag
     end function mxCreateSparse

     logical(c_bool) function mxIsSparse(pm) bind(c, name="mxIsSparse"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsSparse

     type(c_ptr) function mxGetIr(pm) bind(c, name="mxGetIr"//API_VER2)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetIr

     type(c_ptr) function mxGetJc(pm) bind(c, name="mxGetJc"//API_VER2)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetJc

     ! Nonnumeric types
     type(c_ptr) function mxGetData(pm) bind(c, name="mxGetData"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetData

     ! Character
     logical(c_bool) function mxIsChar(pm) bind(c, name="mxIsChar"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsChar

     ! Logical
     logical(c_bool) function mxIsLogical(pm) bind(c, name="mxIsLogical"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsLogical

     logical(c_bool) function mxIsLogicalScalar(array_ptr) bind(c, name="mxIsLogicalScalar"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: array_ptr
     end function mxIsLogicalScalar

     type(c_ptr) function mxCreateLogicalScalar(value) bind(c, name="mxCreateLogicalScalar"//API_VER)
       use iso_c_binding
       import :: mxLogical
       logical(mxLogical), intent(in), value :: value
     end function mxCreateLogicalScalar

     type(c_ptr) function mxGetLogicals_internal(array_ptr) bind(c, name="mxGetLogicals"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: array_ptr
     end function mxGetLogicals_internal

     ! Object
     logical(c_bool) function mxIsClass_internal(pm, classname) bind(c, name="mxIsClass"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
       character(c_char), dimension(*), intent(in) :: classname
     end function mxIsClass_internal

     ! Structure
     type(c_ptr) function mxCreateStructMatrix_internal(m, n, nfields, fieldnames) bind(c, name="mxCreateStructMatrix"//API_VER2)
       use iso_c_binding
       import :: mwSize
       integer(mwSize), intent(in), value :: m, n
       integer(c_int), intent(in), value :: nfields
       type(c_ptr), dimension(*), intent(in) :: fieldnames
     end function mxCreateStructMatrix_internal

     logical(c_bool) function mxIsStruct(pm) bind(c, name="mxIsStruct"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsStruct

     type(c_ptr) function mxGetField_internal(pm, index, fieldname) bind(c, name="mxGetField"//API_VER2)
       use iso_c_binding
       import :: mwIndex
       type(c_ptr), intent(in), value :: pm
       integer(mwIndex), intent(in), value :: index
       character(c_char), dimension(*), intent(in) :: fieldname
     end function mxGetField_internal

     subroutine mxSetField_internal(pm, index, fieldname, pvalue) bind(c, name="mxSetField"//API_VER2)
       use iso_c_binding
       import :: mwIndex
       type(c_ptr), intent(in), value :: pm
       type(c_ptr), intent(in), value :: pvalue
       integer(mwIndex), intent(in), value :: index
       character(c_char), dimension(*), intent(in) :: fieldname
     end subroutine mxSetField_internal

     integer(c_int) function mxGetNumberOfFields(pm) bind(c, name="mxGetNumberOfFields"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxGetNumberOfFields

     ! Cell
     logical(c_bool) function mxIsCell(pm) bind(c, name="mxIsCell"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end function mxIsCell

     type(c_ptr) function mxGetCell_internal(pm, index) bind(c, name="mxGetCell"//API_VER2)
       use iso_c_binding
       import :: mwIndex
       type(c_ptr), intent(in), value :: pm
       integer(mwIndex), intent(in), value :: index
     end function mxGetCell_internal

     !! Delete and Duplicate mxArray
     subroutine mxDestroyArray(pm) bind(c, name="mxDestroyArray"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: pm
     end subroutine mxDestroyArray

     type(c_ptr) function mxDuplicateArray(in) bind(c, name="mxDuplicateArray"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: in
     end function mxDuplicateArray

     !! Convert mxArray

     ! Character
     type(c_ptr) function mxArrayToString_internal(array_ptr) bind(c, name="mxArrayToString"//API_VER)
       use iso_c_binding
       type(c_ptr), intent(in), value :: array_ptr
     end function mxArrayToString_internal
  end interface
contains
  ! Some helper functions to make the interface more Fortran-ish

#if MX_HAS_INTERLEAVED_COMPLEX
  function mxGetDoubles(pm)
    type(c_ptr), intent(in) :: pm
    real(real64), dimension(:), pointer, contiguous :: mxGetDoubles
    call c_f_pointer(mxGetDoubles_internal(pm), mxGetDoubles, [ mxGetNumberOfElements(pm) ])
  end function mxGetDoubles
#endif

  function mxGetPr(pm)
    type(c_ptr), intent(in) :: pm
    real(real64), dimension(:), pointer, contiguous :: mxGetPr
    call c_f_pointer(mxGetPr_internal(pm), mxGetPr, [ mxGetNumberOfElements(pm) ])
  end function mxGetPr

#if MX_HAS_INTERLEAVED_COMPLEX
  function mxGetInt32s(pa)
    type(c_ptr), intent(in) :: pa
    integer(int32), dimension(:), pointer, contiguous :: mxGetInt32s
    call c_f_pointer(mxGetInt32s_internal(pa), mxGetInt32s, [ mxGetNumberOfElements(pa) ])
  end function mxGetInt32s
#endif

#if MX_HAS_INTERLEAVED_COMPLEX
  function mxGetComplexDoubles(pa)
    type(c_ptr), intent(in) :: pa
    complex(real64), dimension(:), pointer, contiguous :: mxGetComplexDoubles
    call c_f_pointer(mxGetComplexDoubles_internal(pa), mxGetComplexDoubles, [ mxGetNumberOfElements(pa) ])
  end function mxGetComplexDoubles
#else
  function mxGetPi(pm)
    type(c_ptr), intent(in) :: pm
    real(real64), dimension(:), pointer, contiguous :: mxGetPi
    call c_f_pointer(mxGetPi_internal(pm), mxGetPi, [ mxGetNumberOfElements(pm) ])
  end function mxGetPi
#endif

  function mxGetLogicals(array_ptr)
    type(c_ptr), intent(in) :: array_ptr
    logical(mxLogical), dimension(:), pointer, contiguous :: mxGetLogicals
    call c_f_pointer(mxGetLogicals_internal(array_ptr), mxGetLogicals, [ mxGetNumberOfElements(array_ptr) ])
  end function mxGetLogicals

  logical(c_bool) function mxIsClass(pm, classname)
    type(c_ptr), intent(in) :: pm
    character(kind=c_char, len=*), intent(in) :: classname
    mxIsclass = mxIsclass_internal(pm, classname // c_null_char)
  end function mxIsClass

  type(c_ptr) function mxCreateStructMatrix(m, n, fieldnames)
    integer(mwSize), intent(in) :: m, n
    character(kind=c_char, len=*), dimension(:), intent(in) :: fieldnames

    character(kind=c_char, len=len(fieldnames)+1), dimension(size(fieldnames)), target :: fieldnames_zero ! Stores zero-terminated strings
    type(c_ptr), dimension(size(fieldnames)) :: fieldnames_ptr ! C arrays of strings
    integer :: i
    do i = 1,size(fieldnames)
       fieldnames_zero(i) = trim(fieldnames(i)) // c_null_char
       fieldnames_ptr(i) = c_loc(fieldnames_zero(i))
    end do
    mxCreateStructMatrix = mxCreateStructMatrix_internal(m, n, size(fieldnames), fieldnames_ptr)
  end function mxCreateStructMatrix

  type(c_ptr) function mxGetField(pm, index, fieldname)
    type(c_ptr), intent(in) :: pm
    integer(mwIndex), intent(in) :: index
    character(kind=c_char, len=*), intent(in) :: fieldname
    mxGetField = mxGetField_internal(pm, index-1, fieldname // c_null_char)
  end function mxGetField

  subroutine mxSetField(pm, index, fieldname, pvalue)
       type(c_ptr), intent(in) :: pm
       type(c_ptr), intent(in) :: pvalue
       integer(mwIndex), intent(in) :: index
       character(kind=c_char, len=*), intent(in) :: fieldname
       call mxSetField_internal(pm, index-1, fieldname // c_null_char, pvalue)
  end subroutine mxSetField

  type(c_ptr) function mxGetCell(pm, index)
    type(c_ptr), intent(in) :: pm
    integer(mwIndex), intent(in) :: index
    mxGetCell = mxGetCell_internal(pm, index-1)
  end function mxGetCell

  function mxArrayToString(pm)
    type(c_ptr), intent(in) :: pm
    character(kind=c_char, len=:), allocatable :: mxArrayToString
    character(kind=c_char), dimension(:), pointer, contiguous :: chararray
    integer :: i
    call c_f_pointer(mxArrayToString_internal(pm), chararray, [ mxGetNumberOfElements(pm) ])
    ! Convert the character array into a character scalar (of length > 1)
    allocate(character(kind=c_char, len=size(chararray)) :: mxArrayToString)
    do i=1,size(chararray)
       mxArrayToString(i:i) = chararray(i)
    end do
  end function mxArrayToString
end module matlab_mat


!!! C MEX API
!!! Listed in same order as https://fr.mathworks.com/help/matlab/call-mex-files-1.html
module matlab_mex
  use matlab_mat
  implicit none

  interface
     integer(c_int) function mexCallMATLAB_internal(nlhs, plhs, nrhs, prhs, functionName) bind(c, name="mexCallMATLAB"//API_VER)
       use iso_c_binding
       integer(c_int), intent(in), value :: nlhs, nrhs
       type(c_ptr), dimension(*), intent(in) :: plhs, prhs
       character(c_char), dimension(*), intent(in) :: functionName
     end function mexCallMATLAB_internal

     type(c_ptr) function mexCallMATLABWithTrap_internal(nlhs, plhs, nrhs, prhs, functionName) &
          bind(c, name="mexCallMATLABWithTrap"//API_VER)
       use iso_c_binding
       integer(c_int), intent(in), value :: nlhs, nrhs
       type(c_ptr), dimension(*), intent(in) :: plhs, prhs
       character(c_char), dimension(*), intent(in) :: functionName
     end function mexCallMATLABWithTrap_internal

     subroutine mexErrMsgTxt_internal(msg) bind(c, name="mexErrMsgTxt"//API_VER)
       use iso_c_binding
       character(c_char), dimension(*), intent(in) :: msg
     end subroutine mexErrMsgTxt_internal

     subroutine mexErrMsgIdAndTxt_internal(id, msg) bind(c, name="mexErrMsgIdAndTxt"//API_VER)
       use iso_c_binding
       character(c_char), dimension(*), intent(in) :: id, msg
     end subroutine mexErrMsgIdAndTxt_internal

     subroutine mexPrintf_internal(message) bind(c, name="mexPrintf"//API_VER)
       use iso_c_binding
       character(c_char), dimension(*), intent(in) :: message
     end subroutine mexPrintf_internal
  end interface
contains
  ! Some helper functions to make the interface more Fortran-ish

  integer(c_int) function mexCallMATLAB(nlhs, plhs, nrhs, prhs, functionName)
    integer(c_int), intent(in) :: nlhs, nrhs
    type(c_ptr), dimension(*), intent(in) :: plhs, prhs
    character(kind=c_char, len=*), intent(in) :: functionName
    mexCallMATLAB = mexCallMATLAB_internal(nlhs, plhs, nrhs, prhs, functionName // c_null_char)
  end function mexCallMATLAB

  type(c_ptr) function mexCallMATLABWithTrap(nlhs, plhs, nrhs, prhs, functionName)
    integer(c_int), intent(in) :: nlhs, nrhs
    type(c_ptr), dimension(*), intent(in) :: plhs, prhs
    character(kind=c_char, len=*), intent(in) :: functionName
    mexCallMATLABWithTrap = mexCallMATLABWithTrap_internal(nlhs, plhs, nrhs, prhs, functionName // c_null_char)
  end function mexCallMATLABWithTrap

  subroutine mexErrMsgTxt(msg)
    character(kind=c_char, len=*), intent(in) :: msg
    call mexErrMsgTxt_internal(msg // c_null_char)
  end subroutine mexErrMsgTxt

  subroutine mexErrMsgIdAndTxt(id, msg)
    character(kind=c_char, len=*), intent(in) :: id, msg
    call mexErrMsgIdAndTxt_internal(id // c_null_char, msg // c_null_char)
  end subroutine mexErrMsgIdAndTxt

  subroutine mexPrintf(message)
    character(kind=c_char, len=*), intent(in) :: message
    call mexPrintf_internal(message // c_null_char)
  end subroutine mexPrintf

  ! Same as mexPrintf(), but trims trailing whitespace, and adds a new line
  subroutine mexPrintf_trim_newline(message)
    character(kind=c_char, len=*), intent(in) :: message
    call mexPrintf_internal(trim(message) // new_line(message) // c_null_char)
  end subroutine mexPrintf_trim_newline
end module matlab_mex
