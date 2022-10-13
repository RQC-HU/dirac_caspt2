import os
import shutil
import pytest
from module_testing import (
    delete_scratch_files,
    is_binary_file_exist,
    create_test_command,
    run_test,
    check_test_returncode,
    get_caspt2_energy_from_output_file,
)


@pytest.mark.dev
def test_c32_co_dev(the_number_of_process: int) -> None:

    # Set file names
    ref_filename = "reference.c32_co_dev.out"  # Reference
    output_filename = "c32_co_dev.caspt2.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.c32_co_dev.caspt2.out"  # latest passed output (After test, the output file is moved to this)

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_file_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    binary_dir = os.path.abspath(os.path.join(test_path, "../../../bin"))  # Set the Built binary directory
    r4dcasci = os.path.abspath(os.path.join(binary_dir, "r4dcascicoexe"))  # CASCI binary
    r4dcaspt2 = os.path.abspath(os.path.join(binary_dir, "r4dcaspt2ocoexe"))  # CASPT2 binary

    # Set delete file list
    delete_files: list[str] = [
        "[A-H]*int*",  # 2-integrals per subspace
        "MDCINTNEW*",  # 2-integrals per MPI process
        "NEWCICOEFF",  # Coefficients of CI
        "CIMAT*",  # CI matrix
        "e0after",  # Energy after CASCI
        "EPS",  # epsilon
        "TRANSFOCK",  # Transformation matrix
        "*mat*",  # Matrix for DMRG
        "fort.*",  # Fortran files
    ]

    # Delete files because of previous test may illegally failed and created files that are not expected
    delete_scratch_files(delete_files, test_path)

    is_binary_file_exist(r4dcasci)
    is_binary_file_exist(r4dcaspt2)

    binaries = [r4dcasci, r4dcaspt2]
    test_command = create_test_command(the_number_of_process, binaries)

    process = run_test(test_command, output_file_path)
    check_test_returncode(process)

    delete_scratch_files(delete_files, test_path)

    ref_energy = get_caspt2_energy_from_output_file(ref_file_path)
    test_energy = get_caspt2_energy_from_output_file(output_file_path)

    # Check whether the output of test run
    # matches the reference to 7th decimal places.
    assert test_energy == pytest.approx(ref_energy, abs=1e-7)

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_file_path)
