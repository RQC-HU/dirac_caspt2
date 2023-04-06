import os
import shutil
import pytest
from module_testing import (
    run_test_caspt2,
    create_test_command_for_ivo,
)


@pytest.mark.dev
def test_ivo_n2_sto3g(mpi_num_process: int, omp_num_threads: int, save: bool) -> None:

    # Set file names
    input_file = "active.ivo.inp"  # Input
    ref_filename = "reference.DFPCMONEW"  # Reference
    test_filename = "DFPCMONEW"  # Test (This file is compared with Reference)
    latest_passed_test = "latest_passed.DFPCMONEW"  # latest passed DFPCMONEW
    output_filename = "c32h_n2_dev.caspt2.out"  # Output
    latest_passed_output = "latest_passed.c32h_n2_dev.caspt2.out"  # latest passed output (After test, the output file is moved to this)

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    test_file_path = os.path.abspath(os.path.join(test_path, test_filename))
    latest_passed_test_file_path = os.path.abspath(os.path.join(test_path, latest_passed_test))
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_output_file_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    binary_dir = os.path.abspath(os.path.join(test_path, "../../../bin"))  # Set the Built binary directory
    dcaspt2 = os.path.join(binary_dir, "dcaspt2")  # Set the dcaspt2 binary path

    test_command = create_test_command_for_ivo(dcaspt2, mpi_num_process, omp_num_threads, input_file, output_file_path, test_path, save)

    run_test_caspt2(test_command)

    # DFPCMONEW format
    # N2                                    Thu Apr   6 11:00:00 2023
    #  2 15 15 54 15 15 54
    #  -0.1075307794799569E+03
    #    -0.0783846159033079    0.0932522354511504   -0.2444662676687903    0.2100050899974797   -0.0207980762625900    0.0061525045631272
    #    -0.0783846159033079    0.0932522354511504   -0.2444662676687903    0.2100050899974797   -0.0207980762625900    0.0061525045631272
    #    -0.0783846159033079    0.0932522354511504   -0.2444662676687903    0.2100050899974797   -0.0207980762625900    0.0061525045631272
    #    -0.0783846159033079    0.0932522354511504   -0.2444662676687903    0.2100050899974797   -0.0207980762625900    0.0061525045631272
    # ...
    #  1 1 1 1 -3 1 1 1 -3 1 1 1 -3 1 1 1 1 1 1 -3 1 1 -3 1 1 1 -3 1 1 1 1 1 1 1 1 -3 1 1 1 -3 1 1 1 -3 1 1 1 1 -3 1 1 -3 1 1 1 -3 1 1 1 1
    #  1 2 3 2

    # Open DFPCMONEW and reference.DFPCMONEW and compare the values (if the values are float, compare the values to 10th decimal places)
    with open(test_file_path, "r") as test_file, open(ref_file_path, "r") as ref_file:
        for test_line, ref_line in zip(test_file, ref_file):
            # if the first value cannot be converted to float, compare the values as strings
            test_values = test_line.split()
            ref_values = ref_line.split()
            try:
                test_float_values = [float(value) for value in test_values]
                ref_float_values = [float(value) for value in ref_values]
                for test_float_value, ref_float_value in zip(test_float_values, ref_float_values):
                    pytest.approx(test_float_value, ref_float_value, abs=1e-10)
            except ValueError:
                assert test_values == ref_values

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_output_file_path)
    shutil.copy(test_file_path, latest_passed_test_file_path)
