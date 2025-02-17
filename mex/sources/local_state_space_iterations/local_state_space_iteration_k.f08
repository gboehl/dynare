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

! Routines and data structures for multithreading over particles in local_state_space_iteration_k
module pparticle
   use iso_c_binding
   use simulation
   use matlab_mex

   implicit none (type, external)

   type tdata
      integer :: nm, nys, endo_nbr, nvar, order, nrestricted, nparticles 
      real(real64), allocatable :: yhat(:,:), e(:,:), ynext(:,:), ys_reordered(:), restrict_var_list(:)
      type(tensor), dimension(:), allocatable :: udr
   end type tdata

   type(tdata) :: thread_data

contains

   subroutine thread_eval(arg) bind(c)
      type(c_ptr), intent(in), value :: arg
      integer, pointer :: im
      integer :: i, j, start, end, q, r, ind
      type(tensor), dimension(:), allocatable :: h
      real(real64), dimension(:), allocatable :: dyu

      ! Checking that the thread number got passed as argument
      if (.not. c_associated(arg)) then
         call mexErrMsgTxt("No argument was passed to thread_eval")
      end if
      call c_f_pointer(arg, im)

      ! Allocating local arrays
      allocate(h(0:thread_data%order), dyu(thread_data%nvar)) 
      do i=0, thread_data%order
         allocate(h(i)%m(thread_data%endo_nbr, thread_data%nvar**i))
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
            thread_data%ynext(i,j) = h(0)%m(ind,1) + thread_data%ys_reordered(ind)
         end do
      end do

   end subroutine thread_eval

end module pparticle

! The code of the local_state_space_iteration_k routine
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
   use pthread
   use pparticle
   implicit none (type, external)

   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs
   type(c_ptr) :: M_mx, options_mx, dr_mx, yhat_mx, epsilon_mx, udr_mx, tmp
   type(tensor), dimension(:), allocatable :: udr
   integer :: order, nstatic, npred, nboth, nfwrd, exo_nbr, endo_nbr, nparticles, nys, nvar, nrestricted, nm
   real(real64), dimension(:), pointer, contiguous :: order_var, ys, restrict_var_list
   real(real64), allocatable :: yhat(:,:), e(:,:), ynext(:,:), ys_reordered(:)
   integer :: i, m, n, rc
   type(c_pthread_t), allocatable :: threads(:)
   integer, allocatable, target :: routines(:)
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
   if (.not. (mxIsDouble(yhat_mx) .and. mxGetM(yhat_mx) >= 1 .and. mxGetN(yhat_mx) >= 1) &
        .or. mxIsComplex(yhat_mx) .or. mxIsSparse(yhat_mx)) then
      call mexErrMsgTxt("1st argument (yhat) should be a real dense vector")
   end if
   if (.not. (mxIsDouble(epsilon_mx) .and. mxGetM(epsilon_mx) >= 1 .or. mxGetN(epsilon_mx) == 1) &
        .or. mxIsComplex(epsilon_mx) .or. mxIsSparse(epsilon_mx)) then
      call mexErrMsgTxt("2nd argument (epsilon) should be a real dense vector")
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
     if (.not. (mxIsDouble(order_var_mx) .and. int(mxGetNumberOfElements(order_var_mx)) == endo_nbr) &
          .or. mxIsComplex(order_var_mx) .or. mxIsSparse(order_var_mx)) then
         call mexErrMsgTxt("Field dr.order_var should be a real dense vector with endo_nbr elements")
      end if
      order_var => mxGetPr(order_var_mx)
   end associate

   associate (ys_mx => mxGetField(dr_mx, 1_mwIndex, "ys"))
      if (.not. (mxIsDouble(ys_mx) .and. int(mxGetNumberOfElements(ys_mx)) == endo_nbr) &
          .or. mxIsComplex(ys_mx) .or. mxIsSparse(ys_mx)) then
         call mexErrMsgTxt("Field dr.ys should be a real dense vector with endo_nbr elements")
      end if
      ys => mxGetPr(ys_mx)
      ! Construct the reordered steady state
      allocate(ys_reordered(endo_nbr))
      do i=1, endo_nbr
         ys_reordered(i) = ys(int(order_var(i)))
      end do
   end associate

   associate (restrict_var_list_mx => mxGetField(dr_mx, 1_mwIndex, "restrict_var_list"))
      if (.not. mxIsDouble(restrict_var_list_mx) .or. mxIsComplex(restrict_var_list_mx) .or. mxIsSparse(restrict_var_list_mx)) then
         call mexErrMsgTxt("Field dr.restrict_var_list should be a real dense vector")
      end if
      nrestricted = size(mxGetPr(restrict_var_list_mx))
      restrict_var_list => mxGetPr(restrict_var_list_mx)
   end associate

   associate (thread_mx => mxGetField(options_mx, 1_mwIndex, "threads"))
      if (.not. c_associated(thread_mx)) then
         call mexErrMsgTxt("Cannot find `threads' in options_")
      end if
      nm = get_int_field(thread_mx, "local_state_space_iteration_k")
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

   allocate(yhat(nys, nparticles), e(exo_nbr, nparticles), ynext(nrestricted, nparticles))
   yhat = reshape(mxGetPr(yhat_mx), [nys, nparticles])
   e = reshape(mxGetPr(epsilon_mx), [exo_nbr, nparticles])


   allocate(udr(0:order)) 
   do i = 0, order
      write (fieldname, '(a2, i1)') "g_", i
      tmp = mxGetField(udr_mx, 1_mwIndex, trim(fieldname))
      if (.not. (c_associated(tmp) .and. mxIsDouble(tmp))) then
         call mexErrMsgTxt(trim(fieldname)//" is not allocated in dr")
      end if
      m = int(mxGetM(tmp))
      n = int(mxGetN(tmp))
      allocate(udr(i)%m(m,n))
      udr(i)%m(1:m,1:n) = reshape(mxGetPr(tmp), [m,n])
   end do

   ! Initializing the global structure containing
   ! useful information for threads
   thread_data%nm = nm
   thread_data%nys = nys
   thread_data%endo_nbr = endo_nbr
   thread_data%nvar = nvar
   thread_data%order = order
   thread_data%nrestricted = nrestricted
   thread_data%nparticles = nparticles
   thread_data%yhat = yhat
   thread_data%e = e
   thread_data%ynext = ynext
   thread_data%udr = udr
   thread_data%ys_reordered = ys_reordered
   thread_data%restrict_var_list = restrict_var_list

   allocate(threads(nm), routines(nm))
   routines = [ (i, i = 1, nm) ]

   if (nm == 1) then
      call thread_eval(c_loc(routines(1)))
   else
      ! Creating the threads
      do i = 1, nm
         rc = c_pthread_create(threads(i), c_null_ptr, c_funloc(thread_eval), c_loc(routines(i)))
      end do

      ! Joining the threads
      do i = 1, nm
         rc = c_pthread_join(threads(i), c_loc(routines(i)))
      end do
   end if

   ! Returning the result
   plhs(1) = mxCreateDoubleMatrix(int(size(restrict_var_list), mwSize), int(nparticles, mwSize), mxREAL)
   mxGetPr(plhs(1)) = reshape(thread_data%ynext, [size(thread_data%ynext)])

end subroutine mexFunction
