! Two buffers allowing the next "panel" of data to arrive
! whilst performing the DGEMM with the existing data
! To perform both send and receive overlapped with compute every iteration,
! one could add a third buffer
!
!mpirun -np 3 cmake_build/summa_overlapping
program summa_1d_overlapping
    use mpi_f08
    use iso_fortran_env, only : dp => real64
    use mpi_utils
    implicit none

    ! --------------------------
    ! Declarations
    ! --------------------------
     ! MPI handles & constants
    integer                      :: rank, np, ierr
    integer,           parameter :: TAG  = 777
    type(MPI_Request)            :: req_recv, req_send

    ! Distribution / indexing metadata
    integer, allocatable :: start_stop_indices(:, :)
    integer, allocatable :: col_sizes(:), displs(:)

    ! Pipelined ring schedule
    integer, allocatable :: recv_src(:), send_dst(:), recv_count(:), owner_recv(:)

    ! Sizes / indices / loop vars
    integer :: ncols_current, offset
    integer :: present_rank, ncols_present
    integer :: j0, j1
    integer :: i, j, irank, s
    logical :: agree, debug_print = .true.

    ! Working buffers & pointer aliases
    real(dp), allocatable, target :: buffer1(:,:), buffer2(:,:)
    real(dp), pointer             :: B_current(:,:), B_next(:,:)
 
    real(dp), parameter :: tol = 1.e-8_dp
    integer,  parameter :: Nrow_a = 4
    integer,  parameter :: Ncol_a = 8
    integer,  parameter :: Nrow_b = 12

    ! Define A on all processes
    real(dp), dimension(Nrow_a, Ncol_a) :: A = reshape([ &
        1._dp,  9._dp, 17._dp, 25._dp,  & ! col 1
        2._dp, 10._dp, 18._dp, 26._dp,  & ! col 2
        3._dp, 11._dp, 19._dp, 27._dp,  & ! col 3
        4._dp, 12._dp, 20._dp, 28._dp,  & ! col 4
        5._dp, 13._dp, 21._dp, 29._dp,  & ! col 5
        6._dp, 14._dp, 22._dp, 30._dp,  & ! col 6
        7._dp, 15._dp, 23._dp, 31._dp,  & ! col 7
        8._dp, 16._dp, 24._dp, 32._dp   & ! col 8
        ], [4,8])

    ! Define a reference for B on all processes, so one does not need to distribute 
    ! with allscatter
    real(dp), dimension(Nrow_b, Ncol_a) :: B_full = reshape([ &
        101._dp, 102._dp, 103._dp, 104._dp, 105._dp, 106._dp, 107._dp, 108._dp, 109._dp, 110._dp, 111._dp, 112._dp, & ! col 1
        201._dp, 202._dp, 203._dp, 204._dp, 205._dp, 206._dp, 207._dp, 208._dp, 209._dp, 210._dp, 211._dp, 212._dp, & ! col 2
        301._dp, 302._dp, 303._dp, 304._dp, 305._dp, 306._dp, 307._dp, 308._dp, 309._dp, 310._dp, 311._dp, 312._dp, & ! col 3
        401._dp, 402._dp, 403._dp, 404._dp, 405._dp, 406._dp, 407._dp, 408._dp, 409._dp, 410._dp, 411._dp, 412._dp, & ! col 4
        501._dp, 502._dp, 503._dp, 504._dp, 505._dp, 506._dp, 507._dp, 508._dp, 509._dp, 510._dp, 511._dp, 512._dp, & ! col 5
        601._dp, 602._dp, 603._dp, 604._dp, 605._dp, 606._dp, 607._dp, 608._dp, 609._dp, 610._dp, 611._dp, 612._dp, & ! col 6
        701._dp, 702._dp, 703._dp, 704._dp, 705._dp, 706._dp, 707._dp, 708._dp, 709._dp, 710._dp, 711._dp, 712._dp, & ! col 7
        801._dp, 802._dp, 803._dp, 804._dp, 805._dp, 806._dp, 807._dp, 808._dp, 809._dp, 810._dp, 811._dp, 812._dp  & ! col 8
        ], [Nrow_b, Ncol_a])

    real(dp), dimension(Nrow_a, Nrow_b) :: C_ref, C
    real(dp), allocatable :: B(:, :)

    ! --------------------------
    ! Main
    ! --------------------------
    call MPI_INIT(ierr)
    call MPI_Comm_rank(mpi_comm_world, rank, ierr)
    call MPI_Comm_size(mpi_comm_world, np,   ierr)

    ! Distribute columns of B
    allocate(start_stop_indices(2, np))
    start_stop_indices = distribute_elements_start_stop(np, size(B_full, 2))

    ! Define B from B_full, such that it is column-distributed
    ncols_current = distribute_elements_from_indices(rank, start_stop_indices)
    allocate(B(Nrow_b, ncols_current))
    offset = start_stop_indices(1, rank + 1)
    do j = 1, ncols_current
        do i = 1, Nrow_b
            B(i, j) = B_full(i, offset + j -1)
        enddo
    enddo

    ! For a single process, perform DGEMM and return
    if (np == 1 .and. ncols_current > 0) then
        call dgemm('N','T', size(A,1), size(C,2), ncols_current, 1.0_dp, &
                A, size(A,1), B, size(B,1), 0.0_dp, C, size(A,1)) 
    end if

    ! All processes need to know the number of columns per process
    allocate(col_sizes(np), displs(np))
    call MPI_Allgather(ncols_current, 1, MPI_INTEGER, col_sizes, 1, MPI_INTEGER, mpi_comm_world, ierr)
    call mpi_displacements(col_sizes, displs)

    ! Buffers and pointer aliases
    allocate(buffer1(Nrow_b, maxval(col_sizes)))
    allocate(buffer2(Nrow_b, maxval(col_sizes)))
    B_current  => buffer1
    B_next => buffer2

    ! Pack local panel into B_current (leading columns)
    ncols_present = ncols_current
    if(ncols_present > 0) B_current(:, 1:ncols_present) = B
    present_rank = rank

    C = 0.0_dp

    ! --- Precompute a branch-free schedule for np iterations ---
    allocate(recv_src(np), send_dst(np), recv_count(np), owner_recv(np))

    ! Owner whose panel will be received at iteration s (1-based)
    owner_recv = [( modulo(rank - s + np, np), s=1, np )]

    ! Source and destination ranks for each iteration 
    recv_src  = [( left_neighbour(rank, np), s=1, np )]
    send_dst  = [( right_neighbour(rank, np), s=1, np )]
    recv_src(np) = MPI_PROC_NULL
    send_dst(np) = MPI_PROC_NULL

    ! Column counts to receive each iteration (last is 0)
    recv_count = [( col_sizes(owner_recv(s)+1), s=1, np )]
    recv_count(np) = 0

    ! --- Fully pipelined ring over np iterations (no per-iter IFs) ---
    do s = 1, np
        ! Post receive for the NEXT panel of this iteration.
        call MPI_Irecv(B_next, size(B_next,1)*recv_count(s), MPI_DOUBLE_PRECISION, &
                        recv_src(s), TAG, MPI_COMM_WORLD, req_recv, ierr)

        ! Compute with the PRESENT panel: C += A(:,j0:j1) * B_current^T
        if (ncols_present > 0) then
            j0 = displs(present_rank+1) + 1
            j1 = j0 + ncols_present - 1
            call dgemm('N','T', size(A,1), size(C,2), ncols_present, 1.0_dp, &
                    A(:,j0:j1), size(A,1), B_current(:,1:ncols_present), size(B_current,1), &
                    1.0_dp, C, size(A,1))
        end if

        ! Send the panel we just used; safe since we’re done reading B_current.
        call MPI_Isend(B_current, size(B_current,1)*ncols_present, MPI_DOUBLE_PRECISION, &
                        send_dst(s), TAG, MPI_COMM_WORLD, req_send, ierr)

        ! Rotate: wait for arrival, swap buffers, update ownership and counts.
        call MPI_Wait(req_recv, MPI_STATUS_IGNORE, ierr)
        call swap_ptr(B_current, B_next)
        present_rank  = owner_recv(s)
        ncols_present = recv_count(s)

        ! Ensure send finished before reusing the buffer next iteration.
        call MPI_Wait(req_send, MPI_STATUS_IGNORE, ierr)
    end do
        
    ! -----------------------------
    ! Validate result
    ! -----------------------------

    ! Reference result
    C_ref = matmul(A, transpose(B_full))

    ! Validate that every process agrees with the reference
    agree = .true.
    do j = 1, size(C, 2)
        do i = 1, size(C, 1)
            agree = abs(C_ref(i, j) - C(i, j)) < tol
        enddo
    enddo

    if (.not. agree) then
        write(*, *) 'C disagrees with the reference on rank ', rank
        error stop 201
    endif

    if (debug_print) then
        do irank = 0, np - 1
            if (rank == irank) then
                write(*, *) 'Distributed C from rank ', irank
                do i = 1, size(C, 1)
                    write(*, *) C(i, :)
                enddo
            endif
            call mpi_barrier(mpi_comm_world, ierr)
        enddo
        if (rank == 0) then
            write(*, *) 'Ref C from rank ', rank
            do i = 1, size(C, 1)
                write(*, *) C_ref(i, :)
            enddo
        endif
    endif

    call MPI_Finalize(ierr)

end program summa_1d_overlapping

