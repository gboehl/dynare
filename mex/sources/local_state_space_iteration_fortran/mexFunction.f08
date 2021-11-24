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

!  input:
!       yhat     values of endogenous variables
!       epsilon  values of the exogenous shock
!       dr       struct containing the folded tensors g_0, g_1, ...
!       M        struct containing the model features
!       options  struct containing the model options
!       udr      struct containing the model unfolded tensors
!  output:
!       ynext    simulated next-period results

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
   use iso_fortran_env
   use iso_c_binding
   use struct
   use matlab_mex
   use partitions
   use simulation
   implicit none

   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs
   type(c_ptr) :: M_mx, options_mx, dr_mx, yhat_mx, epsilon_mx, udr_mx, tmp
   type(pol), dimension(:), allocatable, target :: udr
   integer :: order, nstatic, npred, nboth, nfwrd, exo_nbr, endo_nbr, nparticles, nys, nvar, nrestricted
   real(real64), dimension(:), allocatable :: ys_reordered, dyu
   real(real64), dimension(:), pointer, contiguous :: order_var, ys, restrict_var_list
   real(real64), dimension(:,:), allocatable :: yhat, e, ynext, ynext_all
   type(horner), dimension(:), allocatable :: h
   integer :: i, j, m, n
   character(kind=c_char, len=10) :: fieldname

   yhat_mx = prhs(1)
   epsilon_mx = prhs(2)
   dr_mx = prhs(3)
   M_mx = prhs(4)
   options_mx = prhs(5)
   udr_mx = prhs(6)

   ! Checking the consistence and validity of input arguments
   if (nrhs /= 6 .or. nlhs /= 1) then
      call mexErrMsgTxt("Must have exactly 5 inputs and 1 output")
   end if
   if (.not. (mxIsDouble(yhat_mx) .and. mxGetM(yhat_mx) >= 1 .and. mxGetN(yhat_mx) >= 1)) then
      call mexErrMsgTxt("1st argument (yhat) should be a real vector")
   end if
   if (.not. (mxIsDouble(epsilon_mx) .and. mxGetM(epsilon_mx) >= 1 .or. mxGetN(epsilon_mx) == 1)) then
      call mexErrMsgTxt("2nd argument (epsilon) should be a real vector")
   end if
   if (.not. mxIsStruct(dr_mx)) then
      call mexErrMsgTxt("3rd argument (dr) should be a struct")
   end if
   if (.not. mxIsStruct(M_mx)) then
      call mexErrMsgTxt("4th argument (M) should be a struct")
   end if
   if (.not. mxIsStruct(options_mx)) then
      call mexErrMsgTxt("5th argument (options) should be a struct")
   end if
   if (.not. mxIsStruct(udr_mx)) then
      call mexErrMsgTxt("6th argument (udr) should be a struct")
   end if

   nstatic = get_int_field(M_mx, "nstatic")
   npred = get_int_field(M_mx, "npred")
   nboth = get_int_field(M_mx, "nboth")
   nfwrd = get_int_field(M_mx, "nfwrd")
   endo_nbr = nstatic + npred + nboth + nfwrd
   exo_nbr = get_int_field(M_mx, "exo_nbr")
   order = get_int_field(options_mx, "order")
   nys = npred+nboth
   nvar = nys + exo_nbr

   associate (order_var_mx => mxGetField(dr_mx, 1_mwIndex, "order_var"))
      if (.not. (mxIsDouble(order_var_mx) .and. int(mxGetNumberOfElements(order_var_mx)) == endo_nbr)) then
         call mexErrMsgTxt("Field dr.order_var should be a double precision vector with endo_nbr elements")
      end if
      order_var => mxGetPr(order_var_mx)
   end associate

   associate (ys_mx => mxGetField(dr_mx, 1_mwIndex, "ys"))
      if (.not. (mxIsDouble(ys_mx) .and. int(mxGetNumberOfElements(ys_mx)) == endo_nbr)) then
         call mexErrMsgTxt("Field dr.ys should be a double precision vector with endo_nbr elements")
      end if
      ys => mxGetPr(ys_mx)
      ! Construct the reordered steady state
      allocate(ys_reordered(endo_nbr))
      do i=1, endo_nbr
         ys_reordered(i) = ys(int(order_var(i)))
      end do
   end associate

   associate (restrict_var_list_mx => mxGetField(dr_mx, 1_mwIndex, "restrict_var_list"))
      if (.not. (mxIsDouble(restrict_var_list_mx))) then
         call mexErrMsgTxt("Field dr.restrict_var_list should be a double precision vector")
      end if
      nrestricted = size(mxGetPr(restrict_var_list_mx))
      restrict_var_list => mxGetPr(restrict_var_list_mx)
   end associate

   nparticles = int(mxGetN(yhat_mx));
   if (int(mxGetN(epsilon_mx)) /= nparticles) then
      call mexErrMsgTxt("epsilon and yhat don't have the same number of columns")
   end if
   if (.not. (mxIsDouble(yhat_mx) .and. int(mxGetM(yhat_mx)) == npred + nboth)) then
      call mexErrMsgTxt("yhat should be a double precision matrix with npred+nboth rows")
   end if
   if (.not. (mxIsDouble(epsilon_mx) .and. int(mxGetM(epsilon_mx)) == exo_nbr)) then
      call mexErrMsgTxt("epsilon should be a double precision matrix with exo_nbr rows")
   end if

   allocate(yhat(nys, nparticles), e(exo_nbr, nparticles), ynext(nrestricted, nparticles), ynext_all(endo_nbr, nparticles))
   yhat = reshape(mxGetPr(yhat_mx), [nys, nparticles])
   e = reshape(mxGetPr(epsilon_mx), [exo_nbr, nparticles])

   allocate(h(0:order), udr(0:order)) 
   do i = 0, order
      write (fieldname, '(a2, i1)') "g_", i
      tmp = mxGetField(udr_mx, 1_mwIndex, trim(fieldname))
      if (.not. (c_associated(tmp) .and. mxIsDouble(tmp))) then
         call mexErrMsgTxt(trim(fieldname)//" is not allocated in dr")
      end if
      m = int(mxGetM(tmp))
      n = int(mxGetN(tmp))
      allocate(udr(i)%g(m,n), h(i)%c(endo_nbr, nvar**i))
      udr(i)%g(1:m,1:n) = reshape(mxGetPr(tmp), [m,n])
   end do

   ! Using the Horner algorithm to evaluate the decision rule at the chosen yhat and epsilon
   allocate(dyu(nvar))
   do j=1,nparticles
      dyu(1:nys) = yhat(:,j) 
      dyu(nys+1:) = e(:,j) 
      call eval(h, dyu, udr, endo_nbr, nvar, order)
      ynext_all(:,j) = h(0)%c(:,1) + ys_reordered
      do i=1,nrestricted
         ynext(i,j) = ynext_all(int(restrict_var_list(i)),j) 
      end do
   end do

   plhs(1) = mxCreateDoubleMatrix(int(size(restrict_var_list), mwSize), int(nparticles, mwSize), mxREAL)
   mxGetPr(plhs(1)) = reshape(ynext, [size(ynext)])
 
end subroutine mexFunction
