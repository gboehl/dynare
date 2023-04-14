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
!
!  input:
!       order           the order of approximation, needs order+1 derivatives
!       nstat
!       npred
!       nboth
!       nforw
!       nexog
!       order_moment    order of the moment we need to compute
!       nburn           number of burn-in periods
!       yhat_start      starting value for the full vector of endogenous variables minus its steady-state value
!       shocks          matrix of shocks (number of exogenous variables x number of periods)
!       ysteady         steady-state value for the full vector of endogenous variables
!       dr              structure containing the matrices for the decision rule (g_0, g_1,…)
!  output:
!       mean            estimated `order_moment`-order moment for the full vector of e
!       sim             [optional] matrix of simulations on which the moment was computed

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
   use iso_fortran_env
   use iso_c_binding
   use struct
   use matlab_mex
   use partitions
   use simulation
   implicit none (type, external)

   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs
   type(c_ptr) :: order_mx, nstatic_mx, npred_mx, nboth_mx, nfwrd_mx, nexog_mx, order_moment_mx, nburn_mx, &
   yhat_start_mx, shocks_mx, ysteady_mx, dr_mx, tmp
   type(pol), dimension(:), allocatable, target :: fdr, udr
   integer :: order, nstatic, npred, nboth, nfwrd, exo_nbr, endo_nbr, nys, nvar, nper, nburn, order_moment
   real(real64), dimension(:,:), allocatable :: shocks, sim
   real(real64), dimension(:), allocatable :: dyu, mean
   real(real64), dimension(:), pointer, contiguous :: ysteady, yhat_start
   type(pascal_triangle) :: p
   type(horner), dimension(:), allocatable :: h
   integer :: i, t, d, m, n
   character(kind=c_char, len=10) :: fieldname

   order_mx = prhs(1)
   nstatic_mx = prhs(2)
   npred_mx = prhs(3)
   nboth_mx = prhs(4)
   nfwrd_mx = prhs(5)
   nexog_mx = prhs(6)
   order_moment_mx = prhs(7)
   nburn_mx = prhs(8)
   yhat_start_mx = prhs(9)
   shocks_mx = prhs(10)
   ysteady_mx = prhs(11)
   dr_mx = prhs(12)

   ! Check the number of output arguments
   if (nlhs < 1 .or. nlhs > 2) &
      call mexErrMsgTxt("Must have 1 or 2 outputs")

   ! Checking the consistence and validity of input arguments
   if (nrhs /= 12) then
      call mexErrMsgTxt("Must have exactly 12 inputs")
   end if
   if (.not. (mxIsScalar(order_mx) .and. mxIsNumeric(order_mx))) then
      call mexErrMsgTxt("1st argument (order) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(nstatic_mx) .and. mxIsNumeric(nstatic_mx))) then
      call mexErrMsgTxt("2nd argument (nstat) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(npred_mx) .and. mxIsNumeric(npred_mx))) then
      call mexErrMsgTxt("3rd argument (npred) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(nboth_mx) .and. mxIsNumeric(nboth_mx))) then
      call mexErrMsgTxt("4th argument (nboth) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(nfwrd_mx) .and. mxIsNumeric(nfwrd_mx))) then
      call mexErrMsgTxt("5th argument (nforw) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(nexog_mx) .and. mxIsNumeric(nexog_mx))) then
      call mexErrMsgTxt("6th argument (nexog) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(order_moment_mx) .and. mxIsNumeric(order_moment_mx))) then
      call mexErrMsgTxt("7th argument (order_moment) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(nburn_mx) .and. mxIsNumeric(nburn_mx))) then
      call mexErrMsgTxt("8th argument (nburn) should be a numeric scalar")
   end if
   if (.not. (mxIsDouble(yhat_start_mx) .and. (mxGetM(yhat_start_mx) == 1 .or. mxGetN(yhat_start_mx) == 1)) &
        .or. mxIsComplex(yhat_start_mx) .or. mxIsSparse(yhat_start_mx)) then
      call mexErrMsgTxt("9th argument (yhat_start) should be a real dense vector")
   end if
   if (.not. (mxIsDouble(shocks_mx)) .or. mxIsComplex(shocks_mx) .or. mxIsSparse(shocks_mx)) then
      call mexErrMsgTxt("10th argument (shocks) should be a real dense matrix")
   end if
   if (.not. (mxIsDouble(ysteady_mx) .and. (mxGetM(ysteady_mx) == 1 .or. mxGetN(ysteady_mx) == 1)) &
        .or. mxIsComplex(ysteady_mx) .or. mxIsSparse(ysteady_mx)) then
      call mexErrMsgTxt("11th argument (ysteady) should be a real vector")
   end if
   if (.not. mxIsStruct(dr_mx)) then
      call mexErrMsgTxt("12th argument (dr) should be a struct")
   end if

   ! Converting inputs in Fortran format
   order = int(mxGetScalar(order_mx))
   nstatic = int(mxGetScalar(nstatic_mx))
   npred = int(mxGetScalar(npred_mx))
   nboth = int(mxGetScalar(nboth_mx))
   nfwrd = int(mxGetScalar(nfwrd_mx))
   exo_nbr = int(mxGetScalar(nexog_mx))
   endo_nbr = nstatic+npred+nboth+nfwrd
   nys = npred+nboth
   nvar = nys+exo_nbr
   nburn = int(mxGetScalar(nburn_mx))
   order_moment = int(mxGetScalar(order_moment_mx))

   if (endo_nbr /= int(mxGetM(yhat_start_mx))) then
      call mexErrMsgTxt("yhat_start should have nstat+npred+nboth+nforw rows")
   end if
   yhat_start => mxGetPr(yhat_start_mx)

   if (exo_nbr /= int(mxGetM(shocks_mx))) then
      call mexErrMsgTxt("shocks should have nexog rows")
   end if
   nper = int(mxGetN(shocks_mx))
   allocate(shocks(exo_nbr,nper))
   shocks = reshape(mxGetPr(shocks_mx),[exo_nbr,nper])

   if (.not. (int(mxGetM(ysteady_mx)) == endo_nbr)) then
      call mexErrMsgTxt("ysteady should have nstat+npred+nboth+nforw rows")
   end if
   ysteady => mxGetPr(ysteady_mx)

   allocate(h(0:order), fdr(0:order), udr(0:order)) 
   do i = 0, order
      write (fieldname, '(a2, i1)') "g_", i
      tmp = mxGetField(dr_mx, 1_mwIndex, trim(fieldname))
      if (.not. (c_associated(tmp) .and. mxIsDouble(tmp))) then
         call mexErrMsgTxt(trim(fieldname)//" is not allocated in dr")
      end if
      m = int(mxGetM(tmp))
      n = int(mxGetN(tmp))
      allocate(fdr(i)%g(m,n), udr(i)%g(endo_nbr, nvar**i), h(i)%c(endo_nbr, nvar**i))
      fdr(i)%g(1:m,1:n) = reshape(mxGetPr(tmp), [m,n])
   end do

   udr(0)%g = fdr(0)%g
   udr(1)%g = fdr(1)%g
   if (order > 1) then
      ! Compute the useful binomial coefficients from Pascal's triangle
      p = pascal_triangle(nvar+order-1)
      block
        type(uf_matching), dimension(2:order) :: matching
        ! Pinpointing the corresponding offsets between folded and unfolded tensors
        do d=2,order
           allocate(matching(d)%folded(nvar**d))
           call fill_folded_indices(matching(d)%folded, nvar, d, p)
           udr(d)%g = fdr(d)%g(:,matching(d)%folded)
        end do
      end block
   end if

   allocate(dyu(nvar), mean(endo_nbr), sim(endo_nbr,nper))
   ! Getting the predetermined part of the endogenous variable vector 
   dyu(1:nys) = yhat_start(nstatic+1:nstatic+nys)
   dyu(nys+1:) = shocks(:,1) 
   ! Using the Horner algorithm to evaluate the decision rule at the chosen dyu
   call eval(h, dyu, udr, endo_nbr, nvar, order)
   sim(:,1) = h(0)%c(:,1) + ysteady
   mean = 0.

   ! Carrying out the simulation
   do t=2,nper
      dyu(1:nys) = h(0)%c(nstatic+1:nstatic+nys,1) 
      dyu(nys+1:) = shocks(:,t)
      call eval(h, dyu, udr, endo_nbr, nvar, order)
      sim(:,t) = h(0)%c(:,1) + ysteady
      if (t > nburn) then
         mean = mean + sim(:,t)**order_moment
      end if
   end do
   ! scaling the mean with the number of non-burn-in periods
   mean = mean/(nper-nburn)

   ! Generating output
   plhs(1) = mxCreateDoubleMatrix(int(endo_nbr, mwSize), 1_mwSize, mxREAL)
   mxGetPr(plhs(1)) = mean
   if (nlhs > 1) then
      plhs(2) = mxCreateDoubleMatrix(int(endo_nbr, mwSize), int(nper, mwSize), mxREAL)
      mxGetPr(plhs(2)) = reshape(sim, [size(sim)])
   end if

end subroutine mexFunction
