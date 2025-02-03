program dcaspt2_main
    use module_global_variables, only: rank, nprocs, ierr
    implicit none
    ! MPI initialization and get the number of MPI processes (nprocs) and own process number.
#ifdef HAVE_MPI
    include 'mpif.h'
    real(16)                :: time0, time1
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_rank(MPI_COMM_WORLD, rank, ierr)
    time0 = MPI_Wtime()
#else
    rank = 0; nprocs = 1
#endif
    if (rank == 0) print '(2(A,1X,I0))', 'initialization of mpi, rank :', rank, ' nprocs :', nprocs

    call dcaspt2_run_subprograms

#ifdef HAVE_MPI
    ! MPI finalization
    time1 = MPI_Wtime()
    if (rank == 0) then
        ! Print out the total time (MPI)
        write (*, "(a,e16.6)") "MPI_Wtime :", time1 - time0
    end if
    call MPI_FINALIZE(ierr)
#endif
end program dcaspt2_main
