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

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
   use iso_fortran_env
   use iso_c_binding
   use struct
   use simulation
   use matlab_mex
   use partitions
   implicit none (type, external)

   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs
   type(c_ptr) :: M_mx, options_mx, dr_mx, tmp, g
   type(pol), dimension(:), allocatable, target :: fdr
   integer :: order, npred, nboth, nstatic, nfwrd, endo_nbr, exo_nbr, nys, nvar 
   type(pascal_triangle) :: p
   type(uf_matching), dimension(:), allocatable :: matching 
   character(kind=c_char, len=10), dimension(:), allocatable :: fieldnames
   integer :: d, m, n

   dr_mx = prhs(1)
   M_mx = prhs(2)
   options_mx = prhs(3)

   ! Checking the consistence and validity of input arguments
   if (nrhs /= 3 .or. nlhs /= 1) then
      call mexErrMsgTxt("Must have exactly 3 inputs and 1 output")
   end if
   if (.not. mxIsStruct(dr_mx)) then
      call mexErrMsgTxt("1st argument (dr) should be a struct")
   end if
   if (.not. mxIsStruct(M_mx)) then
      call mexErrMsgTxt("2nd argument (M) should be a struct")
   end if
   if (.not. mxIsStruct(options_mx)) then
      call mexErrMsgTxt("3rd argument (options) should be a struct")
   end if

   nstatic = get_int_field(M_mx, "nstatic")
   npred = get_int_field(M_mx, "npred")
   nboth = get_int_field(M_mx, "nboth")
   nfwrd = get_int_field(M_mx, "nfwrd")
   endo_nbr = nstatic+npred+nboth+nfwrd
   exo_nbr = get_int_field(M_mx, "exo_nbr")
   order = get_int_field(options_mx, "order")
   nys = npred+nboth
   nvar = nys+exo_nbr

   allocate(fdr(0:order), fieldnames(0:order)) 
   do d = 0, order
      write (fieldnames(d), '(a2, i1)') "g_", d
      tmp = mxGetField(dr_mx, 1_mwIndex, trim(fieldnames(d)))
      if (.not. (c_associated(tmp) .and. mxIsDouble(tmp))) then
         call mexErrMsgTxt(trim(fieldnames(d))//" is not allocated in dr")
      end if
      m = int(mxGetM(tmp))
      n = int(mxGetN(tmp))
      allocate(fdr(d)%g(m,n))
      fdr(d)%g = reshape(mxGetPr(tmp), [m,n])
   end do

   plhs(1) = mxCreateStructMatrix(1_mwSize, 1_mwSize, fieldnames)
   g = mxCreateDoubleMatrix(int(endo_nbr, mwSize), 1_mwSize, mxREAL)
   mxGetPr(g) = reshape(fdr(0)%g, [size(fdr(0)%g)])
   call mxSetField(plhs(1), 1_mwIndex, "g_0", g)
   g = mxCreateDoubleMatrix(int(endo_nbr, mwSize), int(nvar, mwSize), mxREAL)
   mxGetPr(g) = reshape(fdr(1)%g, [size(fdr(1)%g)])
   call mxSetField(plhs(1), 1_mwIndex, "g_1", g)

   if (order > 1) then
      ! Compute the useful binomial coefficients from Pascal's triangle
      p = pascal_triangle(nvar+order-1)
      allocate(matching(2:order))
      ! Pinpointing the corresponding offsets between folded and unfolded tensors
      do d=2,order
         allocate(matching(d)%folded(nvar**d))
         call fill_folded_indices(matching(d)%folded, nvar, d, p) 
         g = mxCreateDoubleMatrix(int(endo_nbr, mwSize), int(nvar**d, mwSize), mxREAL)
         mxGetPr(g) = reshape(fdr(d)%g(:,matching(d)%folded), [size(fdr(d)%g(:,matching(d)%folded))])
         call mxSetField(plhs(1), 1_mwIndex, trim(fieldnames(d)), g)
      end do
   end if
 
end subroutine mexFunction
