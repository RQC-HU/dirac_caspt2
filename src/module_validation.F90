module module_validation

    use module_global_variables, only: ndet, rank, nsymrpa, selectroot, nroot
    use module_error, only: stop_with_errorcode
    implicit none
    private
    public validate_ndet

contains
    subroutine validate_ndet()
        implicit none
        if (ndet < selectroot) then
            if (rank == 0) then
                print '(2(a,x,i0,x))', "ERROR: ndet < selectroot :", ndet, "<", selectroot
                print '(a,i0,a)', "Cannot calculate ", selectroot, "th RASCI/CASCI energy"
                print '(a,i0,a,i0)', "because the number of CASCI configuration is ", ndet, " and it is less than ", selectroot
                print *, "Please increase the number of active orbitals or the number of electrons"
                print *, "or decrease the number of selected root."
                print *, "Exit the program."
            end if
            call stop_with_errorcode(1)
        end if

        if (ndet < nroot) then
            if (rank == 0) then
                print '(2(a,x,i0,x))', "WARNING: ndet < nroot:", ndet, "<", nroot
                print '(a,i0,a)', "Cannot print ", nroot, "th RASCI/CASCI energy"
                print '(a,i0,a,i0)', "because the number of CASCI configuration is ", ndet, " and it is less than ", nroot
                print *, "Therefore, replace nroot with the number of CASCI configuration."
                print *, "new nroot = ", ndet
            end if
            nroot = ndet
        end if
    end subroutine validate_ndet
end module module_validation
