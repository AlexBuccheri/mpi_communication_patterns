program summa_1d
    use mpi_f08
    use iso_fortran_env, only : dp => real64
    use mpi_utils
    implicit none

    ! --------------------------
    ! Declarations
    ! --------------------------
    ! MPI
    integer              :: rank, np, ierr, receive_from, send_to
    integer, allocatable :: start_stop_indices(:, :), recvcounts(:), displs(:)
    type(MPI_Request)    :: reqs(2)
    logical              :: receive_from_has_data, this_process_has_data
    integer, parameter :: TAG = 777

    ! Data
    real(dp), parameter :: tol = 1.e-8_dp
    integer,  parameter :: Nrow_a = 4
    integer,  parameter :: Ncol_a = 8
    integer,  parameter :: Nrow_b = 12
    integer             :: Ncol_b_local, offset
 
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
    real(dp), allocatable :: B(:, :), rbuffer(:, :)

    ! Indexing
    integer :: i, j, iring, ireq, irank
    ! Sizes
    integer :: col_start, col_end, ncols_local, size_local
    logical :: agree, debug_print = .true.

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
    ncol_b_local = distribute_elements_from_indices(rank, start_stop_indices)
    allocate(B(Nrow_b, Ncol_b_local))
    offset = start_stop_indices(1, rank + 1)
    do j = 1, Ncol_b_local
        do i = 1, Nrow_b
            B(i, j) = B_full(i, offset + j -1)
        enddo
    enddo

    if (debug_print) then
        do irank = 0, np - 1
            if (rank == irank) then
                write(*, *) 'Distributed B from rank ', irank
                do i = 1, Nrow_b
                    write(*, *) B(i, 1:Ncol_b_local)
                enddo
            endif
            call mpi_barrier(mpi_comm_world, ierr)
        enddo
    endif

    ! All processes need to know the number of columns per process
    allocate(recvcounts(np), displs(np))
    call MPI_Allgather(ncol_b_local, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, mpi_comm_world, ierr)
    call mpi_displacements(recvcounts, displs)

    ! Receive buffer, with upper bound to account for when number of columns has a load
    ! imbalance
    allocate(rbuffer(Nrow_b, maxval(recvcounts)))

    ! ---------------------------
    ! Local contribution
    ! ---------------------------

    ! Guard against if a process has no columns
    if (recvcounts(rank + 1) > 0) then 
        col_start = start_stop_indices(1, rank + 1)
        col_end = start_stop_indices(2, rank + 1)
        ! Initialise C with local data. 
        call dgemm('N','T', size(A, 1), size(C, 2), recvcounts(rank + 1), 1.0_dp, &        
            A(:, col_start:col_end), size(A, 1), B, size(C, 2), 0.0_dp, C, size(A, 1))
    endif
    ! Prevent any other data from touching C until initialisation is done
    call mpi_barrier(mpi_comm_world, ierr)
    
    ! ------------------------------------------------------------
    ! Simplest pairwise rotation: One send/receive per step
    ! Uses a single receive buffer: no rotation, no extra copies but 
    ! no overlapping comms/compute
    ! ------------------------------------------------------------
    do iring = 1, np - 1
        ! Replace these with named functions
        receive_from = modulo(rank - iring + np, np)
        send_to = modulo(rank + iring, np)
        receive_from_has_data = recvcounts(receive_from + 1) > 0
        this_process_has_data = recvcounts(rank + 1) > 0

        ! Receive data from `receive_from` to `rank`
        ireq = 0
        if (receive_from_has_data) then
            ireq = ireq + 1
            ncols_local = recvcounts(receive_from + 1)
            size_local = size(rbuffer, 1) * ncols_local
            call MPI_Irecv(rbuffer(:, 1:ncols_local), size_local, MPI_DOUBLE_PRECISION, receive_from, TAG, mpi_comm_world, reqs(ireq), ierr)
        end if

        ! Send data from rank to `send_to`
        if (this_process_has_data) then
            ireq = ireq + 1
            ncols_local = recvcounts(rank + 1)
            size_local = size(rbuffer, 1) * ncols_local
            call MPI_Isend(B, size_local, MPI_DOUBLE_PRECISION, send_to, TAG, mpi_comm_world, reqs(ireq), ierr)
        end if

        ! No overlapping comms/compute
        if (ireq > 0) call MPI_Waitall(ireq, reqs, MPI_STATUSES_IGNORE, ierr)

        ! Compute contribution from `receive_from`
        if (receive_from_has_data) then
            col_start = start_stop_indices(1, receive_from + 1)
            col_end = start_stop_indices(2, receive_from + 1)
            ncols_local = recvcounts(receive_from + 1)
            call dgemm('N','T', size(A, 1), size(C, 2), ncols_local, 1.0_dp, &
                A(:, col_start:col_end), size(A, 1), rbuffer(:,1:ncols_local), size(C, 2), 1.0_dp, C, size(A, 1))
        end if
    enddo

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
end program
