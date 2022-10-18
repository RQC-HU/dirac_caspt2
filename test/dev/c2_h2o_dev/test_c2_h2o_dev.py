import os
import shutil
import pytest
from module_testing import (
    run_test,
    check_test_returncode,
    create_test_command_for_caspt2,
    get_caspt2_energy_from_output_file,
)


@pytest.mark.dev
def test_c2_h2o_dev(mpi_num_process: int, omp_num_threads: int, save: bool) -> None:

    # Set file names
    input_file = "active.inp"  # Input
    ref_filename = "reference.c2_h2o_dev.out"  # Reference
    output_filename = "c2_h2o_dev.caspt2.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.c2_h2o_dev.caspt2.out"  # latest passed output (After test, the output file is moved to this)

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_file_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    binary_dir = os.path.abspath(os.path.join(test_path, "../../../bin"))  # Set the Built binary directory
    dcaspt2 = os.path.join(binary_dir, "dcaspt2")  # Set the dcaspt2 binary path

    test_command = create_test_command_for_caspt2(dcaspt2, mpi_num_process, omp_num_threads, input_file, output_file_path, test_path, save)
    process = run_test(test_command, output_file_path)
    check_test_returncode(process)

    ref_energy = get_caspt2_energy_from_output_file(ref_file_path)
    test_energy = get_caspt2_energy_from_output_file(output_file_path)

    # Check whether the output of test run
    # matches the reference to 7th decimal places.
    assert test_energy == pytest.approx(ref_energy, abs=1e-8)

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_file_path)
