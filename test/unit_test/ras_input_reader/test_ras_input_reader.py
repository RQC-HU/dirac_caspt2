import os
import shutil
from module_testing import (
    convert_string_list_to_integer_list,
    create_test_command,
    delete_scratch_files,
    get_split_string_list_from_output_file,
    is_binary_file_exist,
    run_test,
)
import pytest


@pytest.mark.dev
def test_ras_input_reader():

    # Set file names
    ref_filename = "expected"  # Reference
    output_filename = "result.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.result.out"  # latest passed output (After test, the output file is moved to this)
    exe_filename = "test_ras_input_reader_exe"  # Executable file

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    is_binary_file_exist(exe_file_path)
    delete_scratch_files([output_filename], test_path)
    test_command = create_test_command(mpi_num_process=1, binaries=[exe_file_path])
    run_test(test_command)

    string_ref = get_split_string_list_from_output_file(ref_file_path)
    ref_int_list = convert_string_list_to_integer_list(string_ref)

    string_result = get_split_string_list_from_output_file(output_file_path)
    result_int_list = convert_string_list_to_integer_list(string_result)

    # Evaluate the difference between references and results
    assert ref_int_list == result_int_list

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)
