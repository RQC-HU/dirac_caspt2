import os
import shutil
from module_testing import (
    create_test_command,
    delete_scratch_files,
    get_stripped_string_from_output_file,
    is_binary_file_exist,
    run_test,
)
import pytest


@pytest.mark.dev
def test_lowercase(env_setup_unittest):
    exe_file_path = env_setup_unittest("test_lowercase_exe")
    print(f"Test command: {exe_file_path}")

    # Set file names
    ref_output_file = "expected"  # Reference
    output_filename = "result.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.result.out"  # latest passed output (After test, the output file is moved to this)

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    ref_output_file_path = os.path.abspath(os.path.join(test_path, ref_output_file))
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))

    is_binary_file_exist(exe_file_path)
    delete_scratch_files([output_filename], test_path)
    test_command = create_test_command(mpi_num_process=1, binaries=[exe_file_path])

    run_test(test_command)

    string_ref = get_stripped_string_from_output_file(ref_output_file_path)
    string_result = get_stripped_string_from_output_file(output_file_path)

    # Evaluate the difference between references and results
    assert string_ref == string_result

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)
