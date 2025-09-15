module mpi_utils
    implicit none
    private

    public :: distribute_elements_from_np, &
              distribute_elements_from_indices, &
              distribute_elements_start_stop, &
              mpi_displacements

contains

  !> @brief Compute the displacements required for (all)gather(v) and (all)scatter(v).
  subroutine mpi_displacements(recvcounts, displs)
    integer, intent(in ) :: recvcounts(:) !< The number of elements that are received from each process
    integer, intent(out) :: displs(:)     !< Entry i specifies the displacement (relative to recvbuf)
    !                                        at which to place the incoming data from process i-1
    integer :: n_processes, i

    n_processes = size(displs)
    if (.not. (size(recvcounts) == n_processes) ) then
        error stop 101
    endif

    displs(1) = 0
    do i = 2, n_processes
      displs(i) = displs(i-1) + recvcounts(i-1)
    enddo

  end subroutine mpi_displacements


   !> @brief Number of elements on process 'rank' for an array of n elements,
   !! distributed over np processes.
   integer function distribute_elements_from_np(np, n, rank)
      integer, intent(in) :: n               !< Number of elements to distribute
      integer, intent(in) :: np              !< Number of processes to distribute across
      integer, intent(in) :: rank            !< Process to compute local number of elements
      integer, allocatable :: indices(:, :)  !< Start and stop indices for each each process
      indices = distribute_elements_start_stop(np, n)
      distribute_elements_from_np = indices(2, rank + 1) - indices(1, rank + 1) + 1
   end function distribute_elements_from_np

   integer function distribute_elements_from_indices(rank, indices)
      integer, intent(in) :: rank            !< Process to compute local number of elements
      integer, intent(in) :: indices(:, :)  !< Start and stop indices for each each process
      distribute_elements_from_indices = indices(2, rank + 1) - indices(1, rank + 1) + 1
   end function distribute_elements_from_indices


   !> @brief Distribute N elements over NP processes.
   !!
   !! Example 1: n = 12 and np = 3
   !! [1, 2, 3, 4 | 5, 6, 7, 8 | 9, 10, 11, 12]
   !!
   !! Example 2: n = 14 and np = 4
   !! If n is not integer-divisible by np, assign the remainders equally between the first n_remainder processes
   !! [1, 2, 3, 4 | 5, 6, 7, 8 |9, 10, 11 |12, 13, 14]
   !!
   function distribute_elements_start_stop(np, n) result(indices)
      integer, intent(in) :: np      !< Number of processes to distribute across
      integer, intent(in) :: n       !< Number of elements to distribute
      integer :: elements_per_process, n_remainder
      integer :: n_assigned_processes, n_assigned_elements
      integer :: iprocess, istart, istop
      integer :: indices(2, np)

      ! When load-balanced
      elements_per_process = int(n / np)
      n_remainder = mod(n, np)

      ! Assign processes with extra elements first
      iprocess = 0
      do while (n_remainder > 0)
         iprocess = iprocess + 1
         istart = 1 + (iprocess - 1)*(elements_per_process + 1)
         istop = iprocess * (elements_per_process + 1)
         indices(:, iprocess) = [istart, istop]
         n_remainder = n_remainder - 1
      end do

      n_assigned_processes = iprocess
      n_assigned_elements = (elements_per_process + 1) * n_assigned_processes

      do iprocess = n_assigned_processes + 1, np
         istart = 1 + ((iprocess - n_assigned_processes - 1) * elements_per_process) + n_assigned_elements
         istop = ((iprocess - n_assigned_processes) * elements_per_process) + n_assigned_elements
         indices(:, iprocess) = [istart, istop]
      end do
 
   end function distribute_elements_start_stop

end module mpi_utils
