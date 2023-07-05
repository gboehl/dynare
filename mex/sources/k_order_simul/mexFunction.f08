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
!       order    the order of approximation, needs order+1 derivatives
!       nstat
!       npred
!       nboth
!       nforw
!       nexog
!       ystart   starting value (full vector of endogenous)
!       shocks   matrix of shocks (nexog x number of period)
!       ysteady  full vector of decision rule's steady
!       dr       structure containing matrices of derivatives (g_0, g_1,…)
!       pruning  boolean stating whether the simulation should be pruned
!  output:
!       res      simulated results

subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
   use iso_c_binding
   use struct
   use simulation
   implicit none (type, external)

   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs
   type(c_ptr) :: order_mx, nstatic_mx, npred_mx, nboth_mx, nfwrd_mx, &
                 &nexog_mx, ystart_mx, shocks_mx, ysteady_mx, dr_mx, &
                 &pruning_mx 
   integer :: order, nstatic, npred, nboth, nfwrd, exo_nbr, endo_nbr, nys, nvar, nper
   real(real64), dimension(:), allocatable :: dy
   real(real64), pointer, contiguous :: ysteady(:), ystart(:), sim(:,:), shocks(:,:)
   logical :: pruning

   order_mx = prhs(1)
   nstatic_mx = prhs(2)
   npred_mx = prhs(3)
   nboth_mx = prhs(4)
   nfwrd_mx = prhs(5)
   nexog_mx = prhs(6)
   ystart_mx = prhs(7)
   shocks_mx = prhs(8)
   ysteady_mx = prhs(9)
   dr_mx = prhs(10)
   pruning_mx = prhs(11)

   ! Checking the consistence and validity of input arguments
   if (nrhs /= 11 .or. nlhs /= 1) then
      call mexErrMsgTxt("Must have exactly 11 inputs and 1 output")
   end if
   if (.not. (mxIsScalar(order_mx)) .and. mxIsNumeric(order_mx)) then
      call mexErrMsgTxt("1st argument (order) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(nstatic_mx)) .and. mxIsNumeric(nstatic_mx)) then
      call mexErrMsgTxt("2nd argument (nstat) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(npred_mx)) .and. mxIsNumeric(npred_mx)) then
      call mexErrMsgTxt("3rd argument (npred) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(nboth_mx)) .and. mxIsNumeric(nboth_mx)) then
      call mexErrMsgTxt("4th argument (nboth) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(nfwrd_mx)) .and. mxIsNumeric(nfwrd_mx)) then
      call mexErrMsgTxt("5th argument (nforw) should be a numeric scalar")
   end if
   if (.not. (mxIsScalar(nexog_mx)) .and. mxIsNumeric(nexog_mx)) then
      call mexErrMsgTxt("6th argument (nexog) should be a numeric scalar")
   end if
   if (.not. (mxIsDouble(ystart_mx) .and. (mxGetM(ystart_mx) == 1 .or. mxGetN(ystart_mx) == 1)) &
        .or. mxIsComplex(ystart_mx) .or. mxIsSparse(ystart_mx)) then
      call mexErrMsgTxt("7th argument (ystart) should be a real dense vector")
   end if
   if (.not. mxIsDouble(shocks_mx) .or. mxIsComplex(shocks_mx) .or. mxIsSparse(shocks_mx)) then
      call mexErrMsgTxt("8th argument (shocks) should be a real dense matrix")
   end if
   if (.not. (mxIsDouble(ysteady_mx) .and. (mxGetM(ysteady_mx) == 1 .or. mxGetN(ysteady_mx) == 1)) &
        .or. mxIsComplex(ysteady_mx) .or. mxIsSparse(ysteady_mx)) then
      call mexErrMsgTxt("9th argument (ysteady) should be a real dense vector")
   end if
   if (.not. mxIsStruct(dr_mx)) then
      call mexErrMsgTxt("10th argument (dr) should be a struct")
   end if
   if (.not. (mxIsLogicalScalar(pruning_mx))) then
      call mexErrMsgTxt("11th argument (pruning) should be a logical scalar")
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
   pruning = mxGetScalar(pruning_mx) == 1._c_double
   nvar = nys+exo_nbr

   if (endo_nbr /= int(mxGetM(ystart_mx))) then
      call mexErrMsgTxt("ystart should have nstat+npred+nboth+nforw rows")
   end if
   ystart => mxGetPr(ystart_mx)

   if (exo_nbr /= int(mxGetM(shocks_mx))) then
      call mexErrMsgTxt("shocks should have nexog rows")
   end if
   nper = int(mxGetN(shocks_mx))
   shocks(1:exo_nbr,1:nper) => mxGetPr(shocks_mx)

   if (.not. (int(mxGetM(ysteady_mx)) == endo_nbr)) then
      call mexErrMsgTxt("ysteady should have nstat+npred+nboth+nforw rows")
   end if
   ysteady => mxGetPr(ysteady_mx)
   ! Initial value for between the states' starting value and the states' 
   ! steady-state value
   dy = ystart(nstatic+1:nstatic+nys)-ysteady(nstatic+1:nstatic+nys)

   if (pruning) then
      dr_mx = mxGetField(dr_mx, 1_mwIndex, "pruning")
      if (.not. mxIsStruct(dr_mx)) then
         call mexErrMsgTxt("dr.pruning should be a struct")
      end if
   end if

   ! Generating output
   plhs(1) = mxCreateDoubleMatrix(int(endo_nbr, mwSize), int(nper+1, mwSize), mxREAL)
   sim(1:endo_nbr,1:(nper+1)) => mxGetPr(plhs(1))
   sim(:,1) = ystart

   if (pruning) then
      call simulate_pruning(sim, dr_mx, ysteady, dy, shocks, order, nstatic, nvar)
   else
      call simulate(sim, dr_mx, ysteady, dy, shocks, order, nstatic, nvar)
   end if

end subroutine mexFunction
