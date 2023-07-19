program main

    use module_file_manager, only: open_formatted_file
    use module_global_variables
    use module_mpi
    use read_input_module, only: read_input
    implicit none

    call mpi_init_wrapper
    if (rank == 0) then
        print '(2(A,1X,I0))', 'initialization of mpi, rank :', rank, ' nprocs :', nprocs
        print *, ''
        print *, ' ENTER DIRAC-CASPT2 PROGRAM written by M. Abe 2007.7.19'
        print *, ''
    end if


    call r4dcasci
    call r4dcaspt2_tra

    call mpi_finalize_wrapper

end program main
