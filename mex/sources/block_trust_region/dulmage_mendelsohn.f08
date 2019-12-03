! Implementation of the Dulmage-Mendelsohn decomposition
!
! Closely follows the description of Pothen and Fan (1990), “Computing the
! Block Triangular Form of a Sparse Matrix”, ACM Transactions on Mathematical
! Software, Vol. 16, No. 4, pp. 303–324
!
! In addition, also computes the directed acyclic graph (DAG) formed by the
! fine blocks, in order to optimize the parallel resolution of a linear system.
! The last (i.e. lower-right) block has no predecessor (since this is an
! *upper* block triangular form), and the other blocks may have predecessors if
! they use variables computed by blocks further on the right.

! Copyright © 2019 Dynare Team
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
! along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

module dulmage_mendelsohn
  use iso_fortran_env
  implicit none

  private
  public :: dm_block, dmperm, dm_blocks

  ! Represents a block in the fine DM decomposition
  type :: dm_block
     integer, dimension(:), allocatable :: row_indices, col_indices
     character :: coarse_type ! Place in the coarse decomposition, either 'H', 'S' or 'V'
     integer, dimension(:), allocatable :: predecessors, successors ! Predecessors and successors in the DAG formed by fine DM blocks
  end type dm_block

  type :: vertex
     integer, dimension(:), allocatable :: edges
     integer :: id ! Index of row of column represented by this vertex
  end type vertex

  type :: matrix_graph
     type(vertex), dimension(:), allocatable :: row_vertices, col_vertices
     integer, dimension(:), allocatable :: row_matching, col_matching ! Maximum matching
     character, dimension(:), allocatable :: row_coarse, col_coarse ! Coarse decomposition
     type(dm_block), dimension(:), allocatable :: fine_blocks ! Fine decomposition
     integer, dimension(:), allocatable :: row_fine_block, col_fine_block ! Maps rows/cols to their block index
   contains
     procedure :: init => matrix_graph_init
     procedure :: maximum_matching => matrix_graph_maximum_matching
     procedure :: coarse_decomposition => matrix_graph_coarse_decomposition
     procedure :: fine_decomposition => matrix_graph_fine_decomposition
     procedure :: fine_blocks_dag => matrix_graph_fine_blocks_dag
  end type matrix_graph

contains
  ! Initialize a matrix_graph object, given a matrix
  subroutine matrix_graph_init(this, mat)
    class(matrix_graph), intent(inout) :: this
    real(real64), dimension(:, :), intent(in) :: mat

    integer :: i, j

    if (size(mat, 1) < size(mat, 2)) &
         error stop "Matrix has less rows than columns; consider transposing it"

    ! Create row vertices
    allocate(this%row_vertices(size(mat, 1)))
    do i = 1, size(mat, 1)
       this%row_vertices(i)%edges = pack([ (j, j=1,size(mat, 2)) ], mat(i, :) /= 0_real64)
       this%row_vertices(i)%id = i
    end do

    ! Create column vertices
    allocate(this%col_vertices(size(mat, 2)))
    do i = 1, size(mat, 2)
       this%col_vertices(i)%edges = pack([ (j, j=1, size(mat, 1)) ], mat(:, i) /= 0_real64)
       this%col_vertices(i)%id = i
    end do
  end subroutine matrix_graph_init

  ! Computes the maximum matching for this graph
  subroutine matrix_graph_maximum_matching(this)
    class(matrix_graph), intent(inout) :: this

    logical, dimension(size(this%row_vertices)) :: row_visited
    integer, dimension(:), allocatable :: col_unmatched, col_unmatched_new
    integer :: i, j

    if (allocated(this%row_matching) .or. allocated(this%col_matching)) &
         error stop "Maximum matching already computed"

    allocate(this%row_matching(size(this%row_vertices)))
    allocate(this%col_matching(size(this%col_vertices)))

    this%row_matching = 0
    this%col_matching = 0
    allocate(col_unmatched(0))

    ! Compute cheap matching
    do i = 1, size(this%col_vertices)
       match_col: do j = 1, size(this%col_vertices(i)%edges)
          associate (v => this%col_vertices(i)%edges(j))
            if (this%row_matching(v) == 0) then
               this%row_matching(v) = i
               this%col_matching(i) = v
               exit match_col
            end if
          end associate
       end do match_col
       if (this%col_matching(i) == 0) col_unmatched = [ col_unmatched, i ]
    end do

    ! Augment matching
    allocate(col_unmatched_new(0))
    do
       row_visited = .false.
       do i = 1, size(col_unmatched)
          call search_augmenting_path(col_unmatched(i))
       end do
       if (size(col_unmatched) == size(col_unmatched_new)) exit
       call move_alloc(from=col_unmatched_new, to=col_unmatched)
       allocate(col_unmatched_new(0))
    end do
  contains
    ! Depth-first search for an augmenting path
    ! The algorithm is written iteratively (and not recursively), because
    ! otherwise the stack could overflow on large matrices
    subroutine search_augmenting_path(col)
      integer, intent(in) :: col

      integer, dimension(:), allocatable :: row_path, col_path ! Path visited so far
      integer, dimension(size(this%col_vertices)) :: col_edge ! For each column, keeps the last visited edge index
      integer :: i

      col_edge = 0
      allocate(col_path(1))
      col_path(1) = col
      allocate(row_path(0))

      do
         ! Depth-first search, starting from last column in the path
         associate (current_col => col_path(size(col_path)))
           col_edge(current_col) = col_edge(current_col) + 1
           if (col_edge(current_col) > size(this%col_vertices(current_col)%edges)) then
              ! We visited all edges from this column, need to backtrack
              if (size(col_path) == 1) then
                 ! The search failed
                 col_unmatched_new = [ col_unmatched_new, col ]
                 return
              end if
              col_path = col_path(:size(col_path)-1)
              row_path = row_path(:size(row_path)-1)
           else
              associate (current_row => this%col_vertices(current_col)%edges(col_edge(current_col)))
                if (.not. row_visited(current_row)) then
                   row_visited(current_row) = .true.
                   row_path = [ row_path, current_row ]
                   if (this%row_matching(current_row) == 0) then
                      ! We found an augmenting path, update matching and return
                      do i = 1, size(col_path)
                         this%row_matching(row_path(i)) = col_path(i)
                         this%col_matching(col_path(i)) = row_path(i)
                      end do
                      return
                   end if
                   col_path = [ col_path, this%row_matching(current_row)]
                end if
              end associate
           end if
         end associate
      end do
    end subroutine search_augmenting_path
  end subroutine matrix_graph_maximum_matching

  ! Computes the coarse decomposition
  ! The {row,col}_coarse arrays are filled with characters 'H', 'S' and 'V'
  subroutine matrix_graph_coarse_decomposition(this)
    class(matrix_graph), intent(inout) :: this

    integer :: i

    if (allocated(this%row_coarse) .or. allocated(this%col_coarse)) &
         error stop "Coarse decomposition already computed"

    if (.not. (allocated(this%row_matching) .and. allocated(this%col_matching))) &
         error stop "Maximum matching not yet computed"

    allocate(this%row_coarse(size(this%row_vertices)))
    allocate(this%col_coarse(size(this%col_vertices)))

    ! Initialize partition
    this%row_coarse = 'S'
    this%col_coarse = 'S'

    ! Identify A_h
    do i = 1, size(this%col_vertices)
       if (this%col_matching(i) /= 0) cycle
       this%col_coarse(i) = 'H'
       call alternating_dfs_col(i)
    end do

    ! Identify A_v
    do i = 1, size(this%row_vertices)
       if (this%row_matching(i) /= 0) cycle
       this%row_coarse(i) = 'V'
       call alternating_dfs_row(i)
    end do
  contains
    ! Depth-first search along alternating path, starting from a column
    ! Written iteratively, to avoid stack overflow on large matrices
    subroutine alternating_dfs_col(col)
      integer, intent(in) :: col

      integer, dimension(size(this%col_vertices)) :: col_edge
      integer, dimension(:), allocatable :: col_path ! path visited so far

      col_edge = 0
      allocate(col_path(1))
      col_path(1) = col

      do
         associate (current_col => col_path(size(col_path)))
           col_edge(current_col) = col_edge(current_col) + 1
           if (col_edge(current_col) > size(this%col_vertices(current_col)%edges)) then
              ! We visited all edges from this column, need to backtrack
              if (size(col_path) == 1) return
              col_path = col_path(:size(col_path)-1)
           else
              associate (current_row => this%col_vertices(current_col)%edges(col_edge(current_col)))
                if (this%row_coarse(current_row) /= 'H') then
                   this%row_coarse(current_row) = 'H'
                   if (this%row_matching(current_row) /= 0) then
                      this%col_coarse(this%row_matching(current_row)) = 'H'
                      col_path = [ col_path, this%row_matching(current_row)]
                   end if
                end if
              end associate
           end if
         end associate
      end do
    end subroutine alternating_dfs_col

    ! Depth-first search along alternating path, starting from a row
    ! This is the perfect symmetric of alternating_dfs_col (swapping "col" and "row")
    subroutine alternating_dfs_row(row)
      integer, intent(in) :: row

      integer, dimension(size(this%row_vertices)) :: row_edge
      integer, dimension(:), allocatable :: row_path ! path visited so far

      row_edge = 0
      allocate(row_path(1))
      row_path(1) = row

      do
         associate (current_row => row_path(size(row_path)))
           row_edge(current_row) = row_edge(current_row) + 1
           if (row_edge(current_row) > size(this%row_vertices(current_row)%edges)) then
              ! We visited all edges from this row, need to backtrack
              if (size(row_path) == 1) return
              row_path = row_path(:size(row_path)-1)
           else
              associate (current_col => this%row_vertices(current_row)%edges(row_edge(current_row)))
                if (this%col_coarse(current_col) /= 'V') then
                   this%col_coarse(current_col) = 'V'
                   if (this%col_matching(current_col) /= 0) then
                      this%row_coarse(this%col_matching(current_col)) = 'V'
                      row_path = [ row_path, this%col_matching(current_col)]
                   end if
                end if
              end associate
           end if
        end associate
      end do
    end subroutine alternating_dfs_row
  end subroutine matrix_graph_coarse_decomposition

  ! Computes the coarse decomposition
  ! Allocate and fills the this%row_fine_* and this%col_fine* arrays
  subroutine matrix_graph_fine_decomposition(this)
    class(matrix_graph), intent(inout) :: this

    ! Index of next block to be discovered
    integer :: block_index

    ! Used in DFS (both in dfs_H_or_V and strongconnect)
    integer, dimension(size(this%row_vertices)) :: row_edge
    integer, dimension(size(this%col_vertices)) :: col_edge

    ! For dfs_H_or_V
    logical, dimension(size(this%row_vertices)) :: row_visited
    logical, dimension(size(this%col_vertices)) :: col_visited

    ! For strongconnect
    integer, dimension(:), allocatable :: col_stack
    logical, dimension(size(this%col_vertices)) :: col_on_stack
    integer, dimension(size(this%col_vertices)) :: col_index, col_lowlink
    integer :: index

    integer :: i

    if (allocated(this%fine_blocks) &
         .or. allocated(this%row_fine_block) .or. allocated(this%col_fine_block)) &
         error stop "Fine decomposition already computed"

    if (.not. (allocated(this%row_coarse) .and. allocated(this%col_coarse))) &
         error stop "Coarse decomposition not yet computed"

    allocate(this%row_fine_block(size(this%row_vertices)))
    allocate(this%col_fine_block(size(this%col_vertices)))
    allocate(this%fine_blocks(0))
    block_index = 1

    row_visited = .false.
    col_visited = .false.
    row_edge = 0
    col_edge = 0

    ! Compute blocks in H
    do i = 1, size(this%row_vertices)
       if (this%row_coarse(i) == 'H' .and. .not. row_visited(i)) call dfs_H_or_V(i, .true., 'H')
    end do
    do i = 1, size(this%col_vertices)
       if (this%col_coarse(i) == 'H' .and. .not. col_visited(i)) call dfs_H_or_V(i, .false., 'H')
    end do

    ! Compute blocks in S (Tarjan algorithm on graph consisting only of column
    ! vertices, where rows have been merged with their matching columns)
    allocate(col_stack(0))
    col_on_stack = .false.
    index = 0
    col_index = -1
    do i = 1, size(this%col_vertices)
       if (this%col_coarse(i) == 'S' .and. col_index(i) == -1) call strongconnect(i)
    end do

    ! Compute blocks in V
    do i = 1, size(this%row_vertices)
       if (this%row_coarse(i) == 'V' .and. .not. row_visited(i)) call dfs_H_or_V(i, .true., 'V')
    end do
    do i = 1, size(this%col_vertices)
       if (this%col_coarse(i) == 'V' .and. .not. col_visited(i)) call dfs_H_or_V(i, .false., 'V')
    end do
  contains
    ! Depth-first search for identifying one block inside the H or V sub-matrix
    subroutine dfs_H_or_V(id, is_row, coarse_type)
      integer, intent(in) :: id            ! Start vertex ID
      logical, intent(in) :: is_row        ! Whether start vertex is a row
      character, intent(in) :: coarse_type ! Either 'H' or 'V'

      integer, dimension(:), allocatable :: path
      type(dm_block) :: new_block

      allocate(path(1))
      path(1) = id

      if (is_row) then
         this%row_fine_block(id) = block_index
         new_block%row_indices = [ id ]
         allocate(new_block%col_indices(0))
         row_visited(id) = .true.
      else
         this%col_fine_block(id) = block_index
         new_block%col_indices = [ id ]
         allocate(new_block%row_indices(0))
         col_visited(id) = .true.
      end if

      do
         associate (current_id => path(size(path)))
           if ((mod(size(path), 2) == 0) .neqv. is_row) then ! Current vertex is a row
              row_edge(current_id) = row_edge(current_id) + 1
              if (row_edge(current_id) > size(this%row_vertices(current_id)%edges)) then
                 if (size(path) == 1) exit
                 path = path(:size(path)-1)
              else
                 associate (next_id => this%row_vertices(current_id)%edges(row_edge(current_id)))
                   if (this%col_coarse(next_id) == coarse_type .and. .not. col_visited(next_id)) then
                      col_visited(next_id) = .true.
                      path = [ path, next_id ]
                      this%col_fine_block(next_id) = block_index
                      new_block%col_indices = [ new_block%col_indices, next_id ]
                   end if
                 end associate
              end if
           else ! Current vertex is a column
              col_edge(current_id) = col_edge(current_id) + 1
              if (col_edge(current_id) > size(this%col_vertices(current_id)%edges)) then
                 if (size(path) == 1) exit
                 path = path(:size(path)-1)
              else
                 associate (next_id => this%col_vertices(current_id)%edges(col_edge(current_id)))
                   if (this%row_coarse(next_id) == coarse_type .and. .not. row_visited(next_id)) then
                      row_visited(next_id) = .true.
                      path = [ path, next_id ]
                      this%row_fine_block(next_id) = block_index
                      new_block%row_indices = [ new_block%row_indices, next_id ]
                   end if
                 end associate
              end if
           end if
         end associate
      end do

      new_block%coarse_type = coarse_type
      this%fine_blocks = [ this%fine_blocks, new_block ]
      block_index = block_index + 1
    end subroutine dfs_H_or_V

    ! Iterative version of the strongconnect routine from:
    ! https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
    subroutine strongconnect(col)
      integer, intent(in) :: col

      integer, dimension(:), allocatable :: col_path

      col_index(col) = index
      col_lowlink(col) = index
      index = index + 1
      col_stack = [ col_stack, col ]
      col_on_stack(col) = .true.

      allocate(col_path(1))
      col_path(1) = col

      do
         associate (current_col => col_path(size(col_path)))
           col_edge(current_col) = col_edge(current_col) + 1
           if (col_edge(current_col) > size(this%col_vertices(current_col)%edges)) then
              ! We visited all edges from this col
              ! If this is a root node, pop the stack and generate an SCC
              if (col_lowlink(current_col) == col_index(current_col)) then
                 associate (current_idx => findloc_back(col_stack, current_col))
                   associate (cols => col_stack(current_idx:))
                     col_on_stack(cols) = .false.
                     associate (rows => this%col_matching(cols))
                       this%col_fine_block(cols) = block_index
                       this%row_fine_block(rows) = block_index
                       this%fine_blocks = [ this%fine_blocks, dm_block(rows, cols, 'S') ]
                     end associate
                   end associate
                   block_index = block_index + 1
                   col_stack = col_stack(:current_idx-1)
                 end associate
              end if

              if (size(col_path) == 1) return

              associate (previous_col => col_path(size(col_path)-1))
                col_lowlink(previous_col) = min(col_lowlink(previous_col), col_lowlink(current_col))
              end associate
              col_path = col_path(:size(col_path)-1)
           else
              associate (next_col => this%row_matching(this%col_vertices(current_col)%edges(col_edge(current_col))))
                if (this%col_coarse(next_col) == 'S' .and. col_index(next_col) == -1) then
                   col_index(next_col) = index
                   col_lowlink(next_col) = index
                   index = index + 1
                   col_stack = [ col_stack, next_col ]
                   col_on_stack(next_col) = .true.

                   col_path = [ col_path, next_col ]
                else if (col_on_stack(next_col)) then
                   col_lowlink(current_col) = min(col_lowlink(current_col), col_index(next_col))
                end if
              end associate
           end if
         end associate
      end do
    end subroutine strongconnect

    ! Equivalent of findloc(array, val, back = .true.) for rank-1 arrays
    ! findloc is in the F2008 standard, but only implemented since gfortran 9
    function findloc_back(array, val)
      integer, dimension(:), intent(in) :: array
      integer, intent(in) :: val
      integer :: findloc_back

      do findloc_back = ubound(array,1),lbound(array,1),-1
         if (array(findloc_back) == val) return
      end do
      error stop "Can’t find element"
    end function findloc_back

  end subroutine matrix_graph_fine_decomposition

  ! Computes the DAG formed by fine blocks
  subroutine matrix_graph_fine_blocks_dag(this)
    class(matrix_graph), intent(inout) :: this

    integer :: blck, i, j
    logical, dimension(size(this%fine_blocks)) :: marked

    ! Compute predecessors in the DAG
    do blck = 1,size(this%fine_blocks)
       marked = .false.
       associate (row_indices => this%fine_blocks(blck)%row_indices)
         do i = 1,size(row_indices)
            associate (row_vertex => this%row_vertices(row_indices(i)))
              do j = 1,size(row_vertex%edges)
                 marked(this%col_fine_block(row_vertex%edges(j))) = .true.
              end do
            end associate
         end do
       end associate
       marked(blck) = .false.
       this%fine_blocks(blck)%predecessors = pack([ (i, i=1,size(this%fine_blocks)) ], marked)
    end do

    ! Compute successors in the DAG
    do blck = 1,size(this%fine_blocks)
       marked = .false.
       associate (col_indices => this%fine_blocks(blck)%col_indices)
         do i = 1,size(col_indices)
            associate (col_vertex => this%col_vertices(col_indices(i)))
              do j = 1,size(col_vertex%edges)
                 marked(this%row_fine_block(col_vertex%edges(j))) = .true.
              end do
            end associate
         end do
       end associate
       marked(blck) = .false.
       this%fine_blocks(blck)%successors = pack([ (i, i=1,size(this%fine_blocks)) ], marked)
    end do
  end subroutine matrix_graph_fine_blocks_dag

  ! Compute the blocks from the fine decomposition, including the DAG they form
  subroutine dm_blocks(mat, blocks)
    real(real64), dimension(:, :), intent(in) :: mat
    type(dm_block), dimension(:), allocatable, intent(out) :: blocks

    type(matrix_graph) :: mg

    call mg%init(mat)
    call mg%maximum_matching
    call mg%coarse_decomposition
    call mg%fine_decomposition
    call mg%fine_blocks_dag

    call move_alloc(mg%fine_blocks, blocks)
  end subroutine dm_blocks

  ! Equivalent of dmperm function of MATLAB/Octave
  subroutine dmperm(mat, row_order, col_order, row_blocks, col_blocks)
    real(real64), dimension(:, :), intent(in) :: mat
    integer, dimension(size(mat, 1)), intent(out) :: row_order
    integer, dimension(size(mat, 2)), intent(out) :: col_order
    integer, dimension(:), allocatable, intent(out) :: row_blocks, col_blocks

    type(matrix_graph) :: mg
    integer :: blck, i

    call mg%init(mat)
    call mg%maximum_matching
    call mg%coarse_decomposition
    call mg%fine_decomposition

    allocate(row_blocks(size(mg%fine_blocks)+1))
    allocate(col_blocks(size(mg%fine_blocks)+1))

    row_blocks(1) = 1
    col_blocks(1) = 1

    do blck = 1,size(mg%fine_blocks)
       associate (row_indices => mg%fine_blocks(blck)%row_indices)
         do i = 1,size(row_indices)
            row_order(row_blocks(blck)+i-1) = row_indices(i)
         end do
         row_blocks(blck+1) = row_blocks(blck)+i-1
       end associate
       associate (col_indices => mg%fine_blocks(blck)%col_indices)
         do i = 1,size(col_indices)
            col_order(col_blocks(blck)+i-1) = col_indices(i)
         end do
         col_blocks(blck+1) = col_blocks(blck)+i-1
       end associate
    end do
  end subroutine dmperm
end module dulmage_mendelsohn
