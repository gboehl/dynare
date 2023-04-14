program dmperm_test
  use iso_fortran_env
  use dulmage_mendelsohn
  implicit none (type, external)

  real(real64), dimension(12, 11) :: M
  integer, dimension(size(M, 1)) :: row_order
  integer, dimension(size(M, 2)) :: col_order
  integer, dimension(:), allocatable :: row_blocks, col_blocks
  type(dm_block), dimension(:), allocatable :: blocks
  integer :: i, j

  ! Matrix given p. 305 of Pothen and Fan (1990)
  M = reshape(real([ &
       1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &
       0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, &
       1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, &
       0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, &
       0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, &
       0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, &
       0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, &
       0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, &
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
       0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, &
       0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, &
       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1  &
       ], real64), shape(M), order = [2, 1])

  ! Shuffle the matrix
  M = M([ 3, 9, 8, 7, 1, 2, 5, 4, 10, 12, 11, 6 ], &
        [ 5, 9, 6, 1, 10, 8, 4, 11, 2, 7, 3 ])

  print *, "Shuffled matrix:"
  do i = 1, size(M, 1)
     write (*, "(11f3.0) ") (M(i, j), j=1,size(M,2))
  end do
  print *

  ! Test dmperm
  call dmperm(M, row_order, col_order, row_blocks, col_blocks)

  print *, "Reordered matrix:"
  do i = 1, size(M, 1)
     write (*, "(11f3.0) ") (M(row_order(i), col_order(j)), j=1,size(M,2))
  end do
  print *

  ! Print fine decomposition
  print *, "Column blocks:", col_blocks
  print *, "Row blocks:", row_blocks
  print *

  ! Test compute_dm_blocks
  call dm_blocks(M, blocks)
  do i = 1, size(blocks)
     print *, "*** Block:", i
     print *, "Coarse type: ", blocks(i)%coarse_type
     print *, "Rows:", blocks(i)%row_indices
     print *, "Columns:", blocks(i)%col_indices
     print *, "Predecessors:", blocks(i)%predecessors
     print *
  end do
end program dmperm_test
