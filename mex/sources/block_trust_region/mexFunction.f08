! Copyright © 2019-2021 Dynare Team
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
  use dulmage_mendelsohn
  use matlab_mex
  use matlab_fcn_closure
  use trust_region
  implicit none

  type(c_ptr), dimension(*), intent(in), target :: prhs
  type(c_ptr), dimension(*), intent(out) :: plhs
  integer(c_int), intent(in), value :: nlhs, nrhs

  real(real64), dimension(:), allocatable, target :: x
  type(dm_block), dimension(:), allocatable, target :: blocks
  integer :: info, i
  real(real64) :: tolf, tolx
  integer :: maxiter
  real(real64), dimension(:), allocatable :: fvec
  real(real64), dimension(:,:), allocatable :: fjac
  logical :: debug, specializedunivariateblocks
  character(len=80) :: debug_msg
  logical(mxLogical), dimension(:), pointer :: isloggedlhs => null(), &
       isauxdiffloggedrhs => null()
  type(c_ptr) :: endo_names, lhs
  logical :: fre ! True if the last block has been solved (i.e. not evaluated), so that residuals must be updated
  integer, dimension(:), allocatable :: evaled_cols ! If fre=.false., lists the columns that have been evaluated so far without updating the residuals

  if (nrhs < 4 .or. nlhs /= 2) then
     call mexErrMsgTxt("Must have at least 7 inputs and exactly 2 outputs")
  end if

  if (.not. ((mxIsChar(prhs(1)) .and. mxGetM(prhs(1)) == 1) .or. mxIsClass(prhs(1), "function_handle"))) then
     call mexErrMsgTxt("First argument (function) should be a string or a function handle")
  end if

  if (.not. (mxIsDouble(prhs(2)) .and. (mxGetM(prhs(2)) == 1 .or. mxGetN(prhs(2)) == 1))) then
     call mexErrMsgTxt("Second argument (initial guess) should be a real vector")
  end if

  if (.not. (mxIsScalar(prhs(3)) .and. mxIsNumeric(prhs(3)))) then
     call mexErrMsgTxt("Third argument (tolf) should be a numeric scalar")
  end if

  if (.not. (mxIsScalar(prhs(4)) .and. mxIsNumeric(prhs(4)))) then
     call mexErrMsgTxt("Fourth argument (tolx) should be a numeric scalar")
  end if

  if (.not. (mxIsScalar(prhs(5)) .and. mxIsNumeric(prhs(5)))) then
     call mexErrMsgTxt("Fifth argument (maxiter) should be a numeric scalar")
  end if

  if (.not. (mxIsLogicalScalar(prhs(6)))) then
     call mexErrMsgTxt("Sixth argument (debug) should be a logical scalar")
  end if

  if (.not. (mxIsStruct(prhs(7)) .and. &
       (mxGetNumberOfFields(prhs(7)) == 0 .or. mxGetNumberOfFields(prhs(7)) == 4))) then
     call mexErrMsgTxt("Seventh argument should be a struct with either 0 or 4 fields")
  end if
  specializedunivariateblocks = (mxGetNumberOfFields(prhs(7)) == 4)

  func => prhs(1)
  tolf = mxGetScalar(prhs(3))
  tolx = mxGetScalar(prhs(4))
  maxiter = int(mxGetScalar(prhs(5)))
  debug = mxGetScalar(prhs(6)) == 1._c_double
  extra_args => prhs(8:nrhs) ! Extra arguments to func are in argument 8 and subsequent ones
  associate (x_mat => mxGetPr(prhs(2)))
    allocate(x(size(x_mat)))
    x = x_mat
  end associate
  if (specializedunivariateblocks) then
     block
       type(c_ptr) :: tmp
       tmp = mxGetField(prhs(7), 1_mwIndex, "isloggedlhs")
       if (.not. (c_associated(tmp) .and. mxIsLogical(tmp) .and. mxGetNumberOfElements(tmp) == size(x))) then
          call mexErrMsgTxt("Seventh argument must have a 'isloggedlhs' field of type logical, of same size as second argument")
       end if
       isloggedlhs => mxGetLogicals(tmp)

       tmp = mxGetField(prhs(7), 1_mwIndex, "isauxdiffloggedrhs")
       if (.not. (c_associated(tmp) .and. mxIsLogical(tmp) .and. mxGetNumberOfElements(tmp) == size(x))) then
          call mexErrMsgTxt("Seventh argument must have a 'isauxdiffloggedrhs' field of type &
               &logical, of same size as second argument")
       end if
       isauxdiffloggedrhs => mxGetLogicals(tmp)

       lhs = mxGetField(prhs(7), 1_mwIndex, "lhs")
       if (.not. (c_associated(lhs) .and. mxIsCell(lhs) .and. mxGetNumberOfElements(lhs) == size(x))) then
          call mexErrMsgTxt("Seventh argument must have a 'lhs' field of type cell, of same size as second argument")
       end if

       endo_names = mxGetField(prhs(7), 1_mwIndex, "endo_names")
       if (.not. (c_associated(endo_names) .and. mxIsCell(endo_names) .and. mxGetNumberOfElements(endo_names) == size(x))) then
          call mexErrMsgTxt("Seventh argument must have a 'endo_names' field of type cell, of same size as second argument")
       end if
     end block

     allocate(evaled_cols(0))
     fre = .false.
  end if

  allocate(fvec(size(x)))
  allocate(fjac(size(x), size(x)))

  ! Compute block decomposition
  nullify(x_indices, f_indices, x_all)
  call matlab_fcn(x, fvec, fjac)
  call dm_blocks(fjac, blocks)

  if (debug) then
     write (debug_msg, "('DYNARE_SOLVE (solve_algo=13|14): number of blocks = ', i0)") size(blocks)
     call mexPrintf_trim_newline(debug_msg)
  end if

  ! Solve the system, starting from bottom-rightmost block
  do i = size(blocks),1,-1
     if (debug) then
        write (debug_msg, "('DYNARE_SOLVE (solve_algo=13|14): solving block ', i0, ' of size ', i0)") &
             i, size(blocks(i)%col_indices)
        call mexPrintf_trim_newline(debug_msg)
     end if

     if (specializedunivariateblocks .and. size(blocks(i)%col_indices) == 1) then
        if (debug) then
           write (debug_msg, "('DYNARE_SOLVE (solve_algo=13|14): solving block ', i0, ' by evaluating RHS')") i
           call mexPrintf_trim_newline(debug_msg)
        end if
        associate (eq => blocks(i)%row_indices(1), var => blocks(i)%col_indices(1))
          if (fre .or. any(abs(fjac(eq, evaled_cols)) > 0._real64)) then
             ! Reevaluation of the residuals is required because the current RHS depends on
             ! variables that potentially have been updated previously.
             nullify(x_indices, f_indices, x_all)
             call matlab_fcn(x, fvec)
             deallocate(evaled_cols) ! This shouldn’t be necessary, but it crashes otherwise with gfortran 8
             allocate(evaled_cols(0))
             fre = .false.
          end if
          evaled_cols = [ evaled_cols, var]
          block
            ! An associate() construct for lhs_eq and endo_name_var makes the
            ! code crash (with double free) using gfortran 8. Hence use a block
            character(kind=c_char, len=:), allocatable :: lhs_eq, endo_name_var
            lhs_eq = mxArrayToString(mxGetCell(lhs, int(eq, mwIndex)))
            endo_name_var = mxArrayToString(mxGetCell(endo_names, int(var, mwIndex)))
            if (lhs_eq == endo_name_var .or. lhs_eq == "log(" // endo_name_var // ")") then
               if (isloggedlhs(eq)) then
                  x(var) = exp(log(x(var)) - fvec(eq))
               else
                  x(var) = x(var) - fvec(eq)
               end if
            else
               if (debug) then
                  write (debug_msg, "('LHS variable is not determined by RHS expression (', i0, ')')") eq
                  call mexPrintf_trim_newline(debug_msg)
                  write (debug_msg, "(a, ' -> ', a)") lhs_eq, endo_name_var
                  call mexPrintf_trim_newline(debug_msg)
               end if
               if (lhs_eq(1:9) == "AUX_DIFF_" .or. lhs_eq(1:13) == "log(AUX_DIFF_") then
                  if (isauxdiffloggedrhs(eq)) then
                     x(var) = exp(log(x(var)) + fvec(eq))
                  else
                     x(var) = x(var) + fvec(eq)
                  end if
               else
                  call mexErrMsgTxt("Algorithm solve_algo=14 cannot be used with this nonlinear problem")
               end if
            end if
          end block
        end associate
        cycle
     else
        if (debug) then
           write (debug_msg, "('DYNARE_SOLVE (solve_algo=13|14): solving block ', i0, ' with trust region routine')") i
           call mexPrintf_trim_newline(debug_msg)
        end if
     end if

     if (size(blocks(i)%col_indices) /= size(blocks(i)%row_indices)) then
        ! Non-square block in DM decomposition
        ! Before erroring out, check whether we are not already at the solution for this block
        ! See also #1851
        if (norm2(fvec(blocks(i)%row_indices)) < tolf) then
           cycle
        else
           call mexErrMsgTxt("DYNARE_SOLVE (solve_algo=13|14): the Dulmage-Mendelsohn &
                &decomposition returned a non-square block. This means that the &
                &Jacobian is singular. You may want to try another value for solve_algo.")
        end if
     end if

     block
       real(real64), dimension(size(blocks(i)%col_indices)) :: x_block
       x_indices => blocks(i)%col_indices
       f_indices => blocks(i)%row_indices
       x_all => x
       x_block = x(x_indices)
       call trust_region_solve(x_block, matlab_fcn, info, tolx, tolf, maxiter)
       x(x_indices) = x_block
     end block

     fre = .true.
  end do

  ! Verify that we have a solution
  ! Note that here we use the ∞-norm, while trust region uses 2-norm; otherwise
  ! this check would almost always fail (because the 2-norm of the full fvec is
  ! larger than the 2-norm of its sub-vectors)
  ! If the check fails, this normally means that the block decomposition was
  ! incorrect (because some element of the Jacobian was numerically zero at the
  ! guess value while not being symbolically zero)
  nullify(x_indices, f_indices, x_all)
  call matlab_fcn(x, fvec)
  if (maxval(abs(fvec)) > tolf) then
     if (debug) &
          call mexPrintf_trim_newline("DYNARE_SOLVE (solve_algo=13|14): residuals still too large, solving for the whole model")
     call trust_region_solve(x, matlab_fcn, info, tolx, tolf, maxiter)
  else
     info = 1
  end if

  plhs(1) = mxCreateDoubleMatrix(int(size(x, 1), mwSize), 1_mwSize, mxREAL)
  mxGetPr(plhs(1)) = x
  if (info == 1) then
     plhs(2) = mxCreateDoubleScalar(0._c_double)
  else
     plhs(2) = mxCreateDoubleScalar(1._c_double)
  end if
end subroutine mexFunction
