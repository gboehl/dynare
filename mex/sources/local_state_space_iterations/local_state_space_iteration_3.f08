! Copyright © 2022-2023 Dynare Team
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

! Routines and data structures for multithreading over particles in local_state_space_iteration_3
module pparticle_3
   use matlab_mex
   use partitions

   implicit none (type, external)

   type tdata_3
      integer :: n, m, s, q, numthreads, xx_size, uu_size, xxx_size, uuu_size
      real(real64), pointer, contiguous :: e(:,:), ghx(:,:), ghu(:,:), &
     &ghxu(:,:), ghxx(:,:), ghuu(:,:), ghs2(:), &
     &ghxxx(:,:), ghuuu(:,:), ghxxu(:,:), ghxuu(:,:), ghxss(:,:), ghuss(:,:), &
     &ss(:), y3(:,:) 
      real(real64), pointer :: yhat3(:,:), yhat2(:,:), yhat1(:,:), ylat3(:,:), &
     &ylat2(:,:), ylat1(:,:)
      type(index), pointer, contiguous :: xx_idcs(:), uu_idcs(:), &
     &xxx_idcs(:), uuu_idcs(:)
      integer, pointer, contiguous :: xx_nbeq(:)
   end type tdata_3

   type(tdata_3) :: td3

contains

   ! Fills y3 as y3 = ybar + ½ghss + ghx·ŷ+ghu·ε + ½ghxx·ŷ⊗ŷ + ½ghuu·ε⊗ε +
   !                  ghxu·ŷ⊗ε + (1/6)·ghxxx ŷ⊗ŷ⊗ŷ + (1/6)·ghuuu·ε⊗ε⊗ε +
   !                  (3/6)·ghxxu·ŷ⊗ŷ⊗ε + (3/6)·ghxuu·ŷ⊗ε⊗ε + 
   !                  (3/6)·ghxss·ŷ + (3/6)·ghuss·ε
   ! in td3
   subroutine thread_eval_3(arg) bind(c)
      type(c_ptr), intent(in), value :: arg
      integer, pointer :: ithread
      integer :: is, im, j, k, start, end, q, r

      ! Checking that the thread number got passed as argument
      if (.not. c_associated(arg)) then
         call mexErrMsgTxt("No argument was passed to thread_eval_3")
      end if
      call c_f_pointer(arg, ithread)

      ! Specifying bounds for the curent thread
      q = td3%s / td3%numthreads
      r = mod(td3%s, td3%numthreads)
      start = (ithread-1)*q+1
      if (ithread < td3%numthreads) then
         end = start+q-1
      else
         end = td3%s
      end if

      do is=start,end
         do im=1,td3%m
            ! y3 = ybar + ½ghss 
            td3%y3(im,is) = td3%ss(im)+0.5_real64*td3%ghs2(im)
            ! y3 += ghx·ŷ+(3/6)·ghxss·ŷ + first n folded indices for ½ghxx·ŷ⊗ŷ
            ! + first n folded indices for (1/6)ghxxx·ŷ⊗ŷ⊗ŷ
            do j=1,td3%n
               td3%y3(im,is) = td3%y3(im,is)+&
              &(0.5_real64*td3%ghxss(j,im)+td3%ghx(j,im))*td3%yhat3(j,is)+&
              &(0.5_real64*td3%ghxx(j,im)+(1._real64/6._real64)*td3%ghxxx(j,im)*td3%yhat3(1, is))*&
              &td3%yhat3(1, is)*td3%yhat3(j,is) 
               ! y3 += ghxu·ŷ⊗ε 
               ! + first n*q folded indices of (3/6)·ghxxu·ŷ⊗ŷ⊗ε
               do k=1,td3%q
                  td3%y3(im,is) = td3%y3(im,is) + &
                 &(td3%ghxu(td3%q*(j-1)+k,im)+&
                 &0.5_real64*td3%ghxxu(td3%q*(j-1)+k,im)*td3%yhat3(1, is))*&
                 &td3%yhat3(j, is)*td3%e(k, is)
               end do
            end do
            ! y3 += ghu·ε+(3/6)·ghuss·ε + first q folded indices of ½ghuu·ε⊗ε
            ! + first q folded indices for (1/6)·ghuuu·ε⊗ε⊗ε
            ! + first n*q folded indices of (3/6)·ghxuu·ŷ⊗ε⊗ε
            do j=1,td3%q
               td3%y3(im,is) = td3%y3(im,is) + &
              &(0.5_real64*td3%ghuss(j,im)+td3%ghu(j,im))*td3%e(j,is) + &
              &(0.5_real64*td3%ghuu(j,im)+(1._real64/6._real64)*td3%ghuuu(j,im)*&
              &td3%e(1, is))*td3%e(1, is)*td3%e(j, is)
               do k=1,td3%n
                  td3%y3(im,is) = td3%y3(im,is) + &
                 &0.5_real64*td3%ghxuu(td3%uu_size*(k-1)+j,im)*&
                 &td3%yhat3(k, is)*td3%e(1, is)*td3%e(j, is)
               end do
            end do
            ! y3 += remaining ½ghxx·ŷ⊗ŷ terms
            ! + the next terms starting from n+1 up to xx_size
            ! of (1/6)ghxxx·ŷ⊗ŷ⊗ŷ
            ! + remaining terms of (3/6)·ghxxu·ŷ⊗ŷ⊗ε
            do j=td3%n+1,td3%xx_size
               td3%y3(im,is) = td3%y3(im,is) + &
              &(0.5_real64*td3%ghxx(j,im)+(1._real64/6._real64)*td3%ghxxx(j,im)*td3%yhat3(1, is))*&
              &td3%yhat3(td3%xx_idcs(j)%coor(1), is)*&
              &td3%yhat3(td3%xx_idcs(j)%coor(2), is)
               do k=1,td3%q
                  td3%y3(im,is) = td3%y3(im,is)+&
                 &0.5_real64*td3%ghxxu(td3%q*(j-1)+k,im)*&
                 &td3%yhat3(td3%xx_idcs(j)%coor(1), is)*&
                 &td3%yhat3(td3%xx_idcs(j)%coor(2), is)*&
                 &td3%e(k, is)
               end do
            end do
            ! y3 += remaining ½ghuu·ε⊗ε terms
            ! + the next uu_size terms starting from q+1
            ! of (1/6)·ghuuu·ε⊗ε⊗ε
            ! + remaining terms of (3/6)·ghxuu·ŷ⊗ε⊗ε
            do j=td3%q+1,td3%uu_size
               td3%y3(im,is) = td3%y3(im,is) + &
              &(0.5_real64*td3%ghuu(j,im)+(1._real64/6._real64)*td3%ghuuu(j,im)*td3%e(1, is))*&
              &td3%e(td3%uu_idcs(j)%coor(1), is)*&
              &td3%e(td3%uu_idcs(j)%coor(2), is)
               do k=1,td3%n
                  td3%y3(im,is) = td3%y3(im,is) + &
                 &0.5_real64*td3%ghxuu(td3%uu_size*(k-1)+j,im)*&
                 &td3%yhat3(k, is)*&
                 &td3%e(td3%uu_idcs(j)%coor(1), is)*&
                 &td3%e(td3%uu_idcs(j)%coor(2), is)
               end do
            end do
            ! y3 += remaining (1/6)·ghxxx·ŷ⊗ŷ⊗ŷ terms
            do j=td3%xx_size+1,td3%xxx_size
               td3%y3(im,is) = td3%y3(im,is)+&
              &(1._real64/6._real64)*td3%ghxxx(j,im)*&
              &td3%yhat3(td3%xxx_idcs(j)%coor(1), is)*&
              &td3%yhat3(td3%xxx_idcs(j)%coor(2), is)*&
              &td3%yhat3(td3%xxx_idcs(j)%coor(3), is)
            end do
            ! y3 += remaining (1/6)ghuuu·ε⊗ε⊗ε terms
            do j=td3%uu_size+1,td3%uuu_size
               td3%y3(im,is) = td3%y3(im,is) + &
              &(1._real64/6._real64)*td3%ghuuu(j,im)*&
              &td3%e(td3%uuu_idcs(j)%coor(1), is)*&
              &td3%e(td3%uuu_idcs(j)%coor(2), is)*&
              &td3%e(td3%uuu_idcs(j)%coor(3), is)
            end do

         end do
      end do

   end subroutine thread_eval_3

   ! Fills ylat1, ylat2, ylat3 and y3 as 
   ! ylat1 = ghx·ŷ1 + ghu·ε 
   ! ylat2 = ½ghss + ghx·ŷ2 + ½ghxx·ŷ1⊗ŷ1 + ½ghuu·ε⊗ε + ghxu·ŷ1⊗ε
   ! ylat3 = ghx·ŷ3 + ghxx·ŷ1⊗ŷ2 + ghxu·ŷ2⊗ε + (1/6)·ghxxx·ŷ1⊗ŷ1⊗ŷ1 
   !         + (1/6)·ghuuu·ε⊗ε⊗ε + (3/6)·ghxxu·ŷ1⊗ŷ1⊗ε 
   !         + (3/6)·ghxuu·ŷ1⊗ε⊗ε + (3/6)·ghxss·ŷ1 + (3/6)·ghuss·ε
   ! y3 = ybar + ylat1 + ylat2 + ylat3
   ! in td3
   subroutine thread_eval_3_pruning(arg) bind(c)
      type(c_ptr), intent(in), value :: arg
      integer, pointer :: ithread
      integer :: is, im, j, k, start, end, q, r, j1, j2
      ! Checking that the thread number got passed as argument
      if (.not. c_associated(arg)) then
         call mexErrMsgTxt("No argument was passed to thread_eval")
      end if
      call c_f_pointer(arg, ithread)

      ! Specifying bounds for the curent thread
      q = td3%s / td3%numthreads
      r = mod(td3%s, td3%numthreads)
      start = (ithread-1)*q+1
      if (ithread < td3%numthreads) then
         end = start+q-1
      else
         end = td3%s
      end if

      do is=start,end
         do im=1,td3%m
            ! y1 = 0 
            ! y2 = ½ghss
            ! y3 = 0
            td3%ylat1(im,is) = td3%ss(im)
            td3%ylat2(im,is) = td3%ss(im)+0.5_real64*td3%ghs2(im)
            td3%ylat3(im,is) = td3%ss(im)
            ! y1 += ghx·ŷ1
            ! y2 += ghx·ŷ2 + first n folded indices for ½ghxx·ŷ1⊗ŷ1
            ! y3 += ghx·ŷ3 +(3/6)·ghxss·ŷ1
            !      + first n folded indices for (1/6)ghxxx·ŷ1⊗ŷ1⊗ŷ1
            do j=1,td3%n
               td3%ylat1(im,is) = td3%ylat1(im,is)+&
              &td3%ghx(j,im)*td3%yhat1(j,is)
               td3%ylat2(im,is) = td3%ylat2(im,is)+&
              &td3%ghx(j,im)*td3%yhat2(j,is)+&
              &0.5_real64*td3%ghxx(j,im)*td3%yhat1(1, is)*td3%yhat1(j, is)
               td3%ylat3(im,is) = td3%ylat3(im,is)+&
              &td3%ghx(j,im)*td3%yhat3(j,is)+&
              &0.5_real64*td3%ghxss(j,im)*td3%yhat1(j,is)+&
              &(1._real64/6._real64)*td3%ghxxx(j,im)*&
              &td3%yhat1(1, is)*td3%yhat1(1, is)*td3%yhat1(j,is) 
               ! y2 += + ghxu·ŷ1⊗ε
               ! y3 += + ghxu·ŷ2⊗ε
               !       + first n*q folded indices of (3/6)·ghxxu·ŷ1⊗ŷ1⊗ε
               do k=1,td3%q
                  td3%ylat2(im,is) = td3%ylat2(im,is)+&
                 &td3%ghxu(td3%q*(j-1)+k,im)*&
                 &td3%yhat1(j, is)*td3%e(k, is)
                  td3%ylat3(im,is) = td3%ylat3(im,is)+&
                 &td3%ghxu(td3%q*(j-1)+k,im)*&
                 &td3%yhat2(j, is)*td3%e(k, is)+&
                 &0.5_real64*td3%ghxxu(td3%q*(j-1)+k,im)*&
                 &td3%yhat1(1, is)*td3%yhat1(j, is)*td3%e(k, is)
               end do
            end do
            ! y1 += + ghu·ε
            ! y2 += + first q folded indices for ½ghuu·ε⊗ε
            ! y3 += + (3/6)·ghuss·ε
            !       + first n*q folded indices of (3/6)·ghxuu·ŷ1⊗ε⊗ε
            !       + first n folded indices of (1/6)·ghuuu·ε⊗ε⊗ε
            do j=1,td3%q
               td3%ylat1(im,is) = td3%ylat1(im,is) + td3%ghu(j,im)*td3%e(j,is)
               td3%ylat2(im,is) = td3%ylat2(im,is) + 0.5_real64*td3%ghuu(j,im)*td3%e(1, is)*td3%e(j, is)
               td3%ylat3(im,is) = td3%ylat3(im,is)+&
              &0.5_real64*td3%ghuss(j,im)*td3%e(j,is)+&
              &(1._real64/6._real64)*td3%ghuuu(j,im)*&
              &td3%e(1, is)*td3%e(1, is)*td3%e(j, is)
               do k=1,td3%n
                  td3%ylat3(im,is) = td3%ylat3(im,is)+&
                 &0.5_real64*td3%ghxuu(td3%uu_size*(k-1)+j,im)*&
                 &td3%yhat1(k, is)*td3%e(1, is)*td3%e(j, is)
               end do
            end do
            ! y2 += remaining ½ghxx·ŷ1⊗ŷ1 terms
            ! y3 += + the next terms starting from n+1 up to xx_size
            !        of (1/6)ghxxx·ŷ1⊗ŷ1⊗ŷ1
            !       + remaining terms of (3/6)·ghxxu·ŷ1⊗ŷ1⊗ε
            do j=td3%n+1,td3%xx_size
               td3%ylat2(im,is) = td3%ylat2(im,is)+&
              &0.5_real64*td3%ghxx(j,im)*&
              &td3%yhat1(td3%xx_idcs(j)%coor(1), is)*&
              &td3%yhat1(td3%xx_idcs(j)%coor(2), is)
               td3%ylat3(im,is) = td3%ylat3(im,is)+&
              &(1._real64/6._real64)*td3%ghxxx(j,im)*&
              &td3%yhat1(1, is)*&
              &td3%yhat1(td3%xx_idcs(j)%coor(1), is)*&
              &td3%yhat1(td3%xx_idcs(j)%coor(2), is)
               do k=1,td3%q
                  td3%ylat3(im,is) = td3%ylat3(im,is)+&
                 &0.5_real64*td3%ghxxu(td3%q*(j-1)+k,im)*&
                 &td3%yhat1(td3%xx_idcs(j)%coor(1), is)*&
                 &td3%yhat1(td3%xx_idcs(j)%coor(2), is)*&
                 &td3%e(k, is)
               end do
            end do
            ! y2 += remaining ½ghuu·ε⊗ε terms
            ! y3 += + remaining terms of (3/6)·ghxuu·ŷ⊗ε⊗ε
            !       + the next uu_size terms starting from q+1
            !         of (1/6)·ghuuu·ε⊗ε⊗ε
            do j=td3%q+1,td3%uu_size
               td3%ylat2(im,is) = td3%ylat2(im,is)+&
              &0.5_real64*td3%ghuu(j,im)*&
              &td3%e(td3%uu_idcs(j)%coor(1), is)*&
              &td3%e(td3%uu_idcs(j)%coor(2), is)
               td3%ylat3(im,is) = td3%ylat3(im,is)+&
              &(1._real64/6._real64)*td3%ghuuu(j,im)*&
              &td3%e(1, is)*&
              &td3%e(td3%uu_idcs(j)%coor(1), is)*&
              &td3%e(td3%uu_idcs(j)%coor(2), is)
               do k=1,td3%n
                  td3%ylat3(im,is) = td3%ylat3(im,is)+&
                 &0.5_real64*td3%ghxuu(td3%uu_size*(k-1)+j,im)*&
                 &td3%yhat1(k, is)*&
                 &td3%e(td3%uu_idcs(j)%coor(1), is)*&
                 &td3%e(td3%uu_idcs(j)%coor(2), is)
               end do
            end do
            ! y3 += remaining (1/6)·ghxxx·ŷ⊗ŷ⊗ŷ terms
            do j=td3%xx_size+1,td3%xxx_size
               td3%ylat3(im,is) = td3%ylat3(im,is)+&
              &(1._real64/6._real64)*td3%ghxxx(j,im)*&
              &td3%yhat1(td3%xxx_idcs(j)%coor(1), is)*&
              &td3%yhat1(td3%xxx_idcs(j)%coor(2), is)*&
              &td3%yhat1(td3%xxx_idcs(j)%coor(3), is)
            end do
            ! y3 += remaining (1/6)ghuuu·ε⊗ε⊗ε terms
            do j=td3%uu_size+1,td3%uuu_size
               td3%ylat3(im,is) = td3%ylat3(im,is)+&
              &(1._real64/6._real64)*td3%ghuuu(j,im)*&
              &td3%e(td3%uuu_idcs(j)%coor(1), is)*&
              &td3%e(td3%uuu_idcs(j)%coor(2), is)*&
              &td3%e(td3%uuu_idcs(j)%coor(3), is)
            end do
            ! y3 += first n folded indices for ghxx·ŷ1⊗ŷ2
            do j=1,td3%xx_size
               j1 = td3%xx_idcs(j)%coor(1)
               j2 = td3%xx_idcs(j)%coor(2)
               if (j1 == j2) then
                  td3%ylat3(im,is) = td3%ylat3(im,is)+&
                 &td3%ghxx(j,im)*&
                 &td3%yhat1(j1,is)*&
                 &td3%yhat2(j2,is)
               else
                  td3%ylat3(im,is) = td3%ylat3(im,is)+&
                 &0.5_real64*td3%ghxx(j,im)*(&
                 &td3%yhat1(j1,is)*&
                 &td3%yhat2(j2,is)+&
                 &td3%yhat1(j2,is)*&
                 &td3%yhat2(j1,is))
               end if
            end do
            td3%y3(im,is) = td3%ylat1(im,is)+td3%ylat2(im,is)+&
                           &td3%ylat3(im,is)-2*td3%ss(im)
         end do
      end do

   end subroutine thread_eval_3_pruning

end module pparticle_3

! The code of the local_state_space_iteration_3 routine
! Input:
! prhs[1] yhat    [double]  n×s array, time t particles if pruning is false
!                           3n×s array, time t particles if pruning is true
! Rows 1 to n contain the pruned first order.
! Rows n+1 to 2*n contain the pruned second order.
! Rows 2*n+1 to 3*n contain the pruned third order.
! prhs[2] e             [double]  q×s array, time t innovations.
! prhs[3] ghx           [double]  m×n array, first order reduced form.
! prhs[4] ghu           [double]  m×q array, first order reduced form.
! prhs[5] ghxx          [double]  m×n² array, second order reduced form.
! prhs[6] ghuu          [double]  m×q² array, second order reduced form.
! prhs[7] ghxu          [double]  m×nq array, second order reduced form.
! prhs[8] ghs2          [double]  m×1 array, second order reduced form.
! prhs[9] ghxxx         [double]  m×n array, third order reduced form.
! prhs[10] ghuuu        [double]  m×q array, third order reduced form.
! prhs[11] ghxxu        [double]  m×n²q array, third order reduced form.
! prhs[12] ghxuu        [double]  m×nq² array, third order reduced form.
! prhs[13] ghxss        [double]  m×n array, third order reduced form.
! prhs[14] ghuss        [double]  m×q array, third order reduced form.
! prhs[15] ss           [double]  m×1 array, deterministic steady state
! prhs[16] numthreads   [double]  num of threads
! prhs[17] pruning      [double]  pruning option
!
! Output:
! plhs[1] y3             [double]  m×s array, time t+1 particles.
! plhs[2] ylat           [double]  3m×s array, time t+1 particles for the 
! pruning latent variables up to the 3rd order. Rows 1 to m contain the pruned
! first order. Rows m+1 to 2*m contain the pruned second order. Rows 2*m+1
! to 3*m contain the pruned third order.
subroutine mexFunction(nlhs, plhs, nrhs, prhs) bind(c, name='mexFunction')
   use iso_c_binding
   use matlab_mex
   use pascal
   use partitions
   use pthread
   use pparticle_3
   implicit none (type, external)

   type(c_ptr), dimension(*), intent(in), target :: prhs
   type(c_ptr), dimension(*), intent(out) :: plhs
   integer(c_int), intent(in), value :: nlhs, nrhs
   integer :: n, m, s, q, numthreads
   real(real64), pointer, contiguous :: ghx(:,:), ghu(:,:), ghxx(:,:), &
  &ghuu(:,:), ghxu(:,:), ghxxx(:,:), ghuuu(:,:), ghxxu(:,:),  &
  &ghxuu(:,:), ghxss(:,:), ghuss(:,:), yhatlat(:,:), ylat(:,:)
   integer :: i, j, k, xx_size, uu_size, xxx_size, uuu_size, rc
   character(kind=c_char, len=10) :: arg_nber
   type(pascal_triangle) :: p
   integer, allocatable :: xxx_nbeq(:), &
  &uu_nbeq(:), uuu_nbeq(:), xx_off(:), uu_off(:), &
  &xxx_off(:), uuu_off(:)
   integer, allocatable, target :: xx_nbeq(:)
   type(c_pthread_t), allocatable :: threads(:)
   integer, allocatable, target :: routines(:)
   logical :: pruning

   ! 0. Checking the consistency and validity of input arguments
   if (nrhs /= 17) then
      call mexErrMsgTxt("Must have exactly 17 inputs")
   end if

   if (nlhs > 2) then
      call mexErrMsgTxt("Too many output arguments.")
   end if

   do i=1,15
      if (.not. (c_associated(prhs(i)) .and. mxIsDouble(prhs(i)) .and. &
          (.not. mxIsComplex(prhs(i))) .and. (.not. mxIsSparse(prhs(i))))) then
         write (arg_nber,"(i2)") i 
         call mexErrMsgTxt("Argument " // trim(arg_nber) // " should be a real dense matrix")
      end if
   end do

   if (.not. (c_associated(prhs(16)) .and. mxIsScalar(prhs(16)) .and. &
       mxIsNumeric(prhs(16)))) then
      call mexErrMsgTxt("Argument 16 should be a numeric scalar")
   end if
   numthreads = int(mxGetScalar(prhs(16)))
   if (numthreads <= 0) then
      call mexErrMsgTxt("Argument 16 should be a positive integer")
   end if
   td3%numthreads = numthreads

   if (.not. (c_associated(prhs(17)) .and. mxIsLogicalScalar(prhs(17)))) then
      call mexErrMsgTxt("Argument 17 should be a logical scalar")
   end if
   pruning = mxGetScalar(prhs(17)) == 1._c_double

   if (pruning) then
      n = int(mxGetM(prhs(1)))/3   ! Number of states.
   else
      n = int(mxGetM(prhs(1)))     ! Number of states.
   end if
   s = int(mxGetN(prhs(1)))   ! Number of particles.
   q = int(mxGetM(prhs(2)))   ! Number of innovations.
   m = int(mxGetM(prhs(3)))   ! Number of elements in the union of states and observed variables.
   td3%n = n
   td3%s = s
   td3%q = q
   td3%m = m

   if ((s /= mxGetN(prhs(2)))            &  ! Number of columns for epsilon
      &.or. (n /= mxGetN(prhs(3)))       &  ! Number of columns for ghx
      &.or. (m /= mxGetM(prhs(4)))       &  ! Number of rows for ghu
      &.or. (q /= mxGetN(prhs(4)))       &  ! Number of columns for ghu
      &.or. (m /= mxGetM(prhs(5)))       &  ! Number of rows for ghxx
      &.or. (n*n /= mxGetN(prhs(5)))     &  ! Number of columns for ghxx
      &.or. (m /= mxGetM(prhs(6)))       &  ! Number of rows for ghuu
      &.or. (q*q /= mxGetN(prhs(6)))     &  ! Number of columns for ghuu
      &.or. (m /= mxGetM(prhs(7)))       &  ! Number of rows for ghxu
      &.or. (n*q /= mxGetN(prhs(7)))     &  ! Number of columns for ghxu
      &.or. (m /= mxGetM(prhs(8)))       &  ! Number of rows for ghs2
      &.or. (m /= mxGetM(prhs(9)))       &  ! Number of rows for ghxxx
      &.or. (n*n*n /= mxGetN(prhs(9)))   &  ! Number of columns for ghxxx
      &.or. (m /= mxGetM(prhs(10)))      &  ! Number of rows for ghuuu
      &.or. (q*q*q /= mxGetN(prhs(10)))  &  ! Number of columns for ghuuu
      &.or. (m /= mxGetM(prhs(11)))      &  ! Number of rows for ghxxu
      &.or. (n*n*q /= mxGetN(prhs(11)))  &  ! Number of columns for ghxxu
      &.or. (m /= mxGetM(prhs(12)))      &  ! Number of rows for ghxuu
      &.or. (n*q*q /= mxGetN(prhs(12)))  &  ! Number of columns for ghxuu
      &.or. (m /= mxGetM(prhs(13)))      &  ! Number of rows for ghxss
      &.or. (n /= mxGetN(prhs(13)))      &  ! Number of columns for ghxss
      &.or. (m /= mxGetM(prhs(14)))      &  ! Number of rows for ghuss
      &.or. (q /= mxGetN(prhs(14)))      &  ! Number of columns for ghuss
      &.or. (m /= mxGetM(prhs(15)))      &  ! Number of rows for ss
      &) then
      call mexErrMsgTxt("Input dimension mismatch")
   end if

   ! 1. Getting relevant information to take advantage of symmetries
   ! There are symmetries in the ghxx, ghuu, ghxxx, ghuuu, ghxxu and ghxuu terms
   ! that we may exploit to avoid unnecessarily repeating operations in matrix-vector
   ! multiplications, e.g in ghxx·ŷ⊗ŷ. 
   ! In matrix-vector multiplications such as ghxx·ŷ⊗ŷ, we loop through all the folded offsets
   ! and thus need for each one of them :
   !    (i) the corresponding folded index, e.g (α₁,α₂), α₁≤α₂ for ghxx
   !    (i) the corresponding offset in the unfolded matrix
   !    (ii) the corresponding number of equivalent unfolded indices (1 if α₁=α₂, 2 otherwise)
   ! It is better to compute these beforehand as it avoids repeating the calculation for 
   ! each particle. The `folded_offset_loop` routine carries out this operation. 

   p = pascal_triangle(max(n,q)+3-1)
   xx_size = get(2,n+2-1,p)
   uu_size = get(2,q+2-1,p)
   xxx_size = get(3,n+3-1,p)
   uuu_size = get(3,q+3-1,p)

   td3%xx_size = xx_size
   td3%uu_size = uu_size
   td3%xxx_size = xxx_size
   td3%uuu_size = uuu_size

   allocate(td3%xx_idcs(xx_size), td3%uu_idcs(uu_size), &
           &td3%xxx_idcs(xxx_size), td3%uuu_idcs(uuu_size), &
           &xx_off(xx_size), uu_off(uu_size), &
           &xxx_off(xxx_size), uuu_off(uuu_size), &
           &xx_nbeq(xx_size), uu_nbeq(uu_size), &
           &xxx_nbeq(xxx_size), uuu_nbeq(uuu_size))

   call folded_offset_loop(td3%xx_idcs, xx_nbeq, &
                          &xx_off, n, 2, p)
   call folded_offset_loop(td3%uu_idcs, uu_nbeq, &
                          &uu_off, q, 2, p)
   call folded_offset_loop(td3%xxx_idcs, xxx_nbeq, &
                          &xxx_off, n, 3, p)
   call folded_offset_loop(td3%uuu_idcs, uuu_nbeq, &
                          &uuu_off, q, 3, p)

   ! 1. Storing the relevant input variables in Fortran
   if (pruning) then
      yhatlat(1:3*n,1:s) => mxGetPr(prhs(1))
      td3%yhat1 => yhatlat(1:n,1:s)
      td3%yhat2 => yhatlat(n+1:2*n,1:s)
      td3%yhat3 => yhatlat(2*n+1:3*n,1:s)
      ! td3%xx_nbeq => xx_nbeq
   else
      td3%yhat3(1:n,1:s) => mxGetPr(prhs(1))
   end if
   td3%e(1:q,1:s) => mxGetPr(prhs(2))
   ghx(1:m,1:n) => mxGetPr(prhs(3))
   ghu(1:m,1:q) => mxGetPr(prhs(4))
   ghxx(1:m,1:n*n) => mxGetPr(prhs(5))
   ghuu(1:m,1:q*q) => mxGetPr(prhs(6))
   ghxu(1:m,1:n*q) => mxGetPr(prhs(7))
   td3%ghs2 => mxGetPr(prhs(8))
   ghxxx(1:m,1:n*n*n) => mxGetPr(prhs(9))
   ghuuu(1:m,1:q*q*q) => mxGetPr(prhs(10))
   ghxxu(1:m,1:n*n*q) => mxGetPr(prhs(11))
   ghxuu(1:m,1:n*q*q) => mxGetPr(prhs(12))
   ghxss(1:m,1:n) => mxGetPr(prhs(13))
   ghuss(1:m,1:q) => mxGetPr(prhs(14))
   td3%ss => mxGetPr(prhs(15))

   ! Getting a transposed folded copy of the unfolded tensors
   ! for future loops to be more efficient
   allocate(td3%ghx(n,m), td3%ghu(q,m),&
           &td3%ghuu(uu_size,m), td3%ghxu(n*q,m), &
           &td3%ghxx(xx_size,m), &
           &td3%ghxxx(xxx_size,m), td3%ghuuu(uuu_size,m), &
           &td3%ghxxu(xx_size*q,m), td3%ghxuu(n*uu_size,m), &
           &td3%ghxss(n,m), td3%ghuss(q,m))
   do i=1,m
      do j=1,n
         td3%ghx(j,i) = ghx(i,j)
         td3%ghxss(j,i) = ghxss(i,j)
         td3%ghxx(j,i) = xx_nbeq(j)*ghxx(i,xx_off(j))
         td3%ghxxx(j,i) = xxx_nbeq(j)*ghxxx(i,xxx_off(j))
         do k=1,q
            td3%ghxu(q*(j-1)+k,i) = ghxu(i,q*(j-1)+k) 
            td3%ghxxu(q*(j-1)+k,i) = xx_nbeq(j)*ghxxu(i,q*(xx_off(j)-1)+k)
         end do
      end do
      do j=n+1,xx_size
         td3%ghxx(j,i) = xx_nbeq(j)*ghxx(i,xx_off(j))
         td3%ghxxx(j,i) = xxx_nbeq(j)*ghxxx(i,xxx_off(j))
         do k=1,q
            td3%ghxxu(q*(j-1)+k,i) = xx_nbeq(j)*ghxxu(i,q*(xx_off(j)-1)+k)
         end do
      end do
      do j=xx_size+1,xxx_size
         td3%ghxxx(j,i) = xxx_nbeq(j)*ghxxx(i,xxx_off(j))
      end do
      do j=1,q
         td3%ghu(j,i) = ghu(i,j)
         td3%ghuss(j,i) = ghuss(i,j)
         td3%ghuu(j,i) = uu_nbeq(j)*ghuu(i,uu_off(j))
         td3%ghuuu(j,i) = uuu_nbeq(j)*ghuuu(i,uuu_off(j))
         do k=1,n
            td3%ghxuu(uu_size*(k-1)+j,i) = uu_nbeq(j)*ghxuu(i,q*q*(k-1)+uu_off(j))
         end do
      end do
      do j=q+1,uu_size
         td3%ghuu(j,i) = uu_nbeq(j)*ghuu(i,uu_off(j))
         td3%ghuuu(j,i) = uuu_nbeq(j)*ghuuu(i,uuu_off(j))
         do k=1,n
            td3%ghxuu(uu_size*(k-1)+j,i) = uu_nbeq(j)*ghxuu(i,q*q*(k-1)+uu_off(j))
         end do
      end do
      do j=uu_size+1,uuu_size
         td3%ghuuu(j,i) = uuu_nbeq(j)*ghuuu(i,uuu_off(j))
      end do
   end do

   ! 3. Implementing the calculations:

   plhs(1) = mxCreateDoubleMatrix(int(m, mwSize), int(s, mwSize), mxREAL)
   td3%y3(1:m,1:s) => mxGetPr(plhs(1))
   if (pruning) then
      plhs(2) = mxCreateDoubleMatrix(int(3*m, mwSize), int(s, mwSize), mxREAL)
      ylat(1:3*m,1:s) => mxGetPr(plhs(2))
      td3%ylat1 => ylat(1:m,1:s)
      td3%ylat2 => ylat(m+1:2*m,1:s)
      td3%ylat3 => ylat(2*m+1:3*m,1:s)
   end if

   allocate(threads(numthreads), routines(numthreads))
   routines = [ (i, i = 1, numthreads) ]

   if (numthreads == 1) then
      if (pruning) then
         call thread_eval_3_pruning(c_loc(routines(1)))
      else 
         call thread_eval_3(c_loc(routines(1)))
      end if
   else
      ! Creating the threads
      if (pruning) then
         do i = 1, numthreads
            rc = c_pthread_create(threads(i), c_null_ptr, c_funloc(thread_eval_3_pruning), c_loc(routines(i)))
         end do
      else
         do i = 1, numthreads
            rc = c_pthread_create(threads(i), c_null_ptr, c_funloc(thread_eval_3), c_loc(routines(i)))
         end do
      end if

      ! Joining the threads
      do i = 1, numthreads
         rc = c_pthread_join(threads(i), c_loc(routines(i)))
      end do
   end if

end subroutine mexFunction
