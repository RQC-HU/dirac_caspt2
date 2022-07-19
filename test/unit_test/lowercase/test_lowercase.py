import os
import shutil
from module_testing import (
    is_binary_file_exist,
    create_test_command,
    run_test,
    check_test_returncode,
    get_stripped_string_from_output_file,
)


def test_lowercase():

    # Set file names
    ref_filename = "expected"  # Reference
    output_filename = "result.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.result.out"  # latest passed output (After test, the output file is moved to this)
    exe_filename = "test_lowercase_exe"  # Executable file

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    test_command = create_test_command(the_number_of_process=1, binaries=[exe_file_path])

    process = run_test(test_command, output_file_path)
    check_test_returncode(process)

    string_ref = get_stripped_string_from_output_file(ref_file_path)
    string_result = get_stripped_string_from_output_file(output_file_path)

    # Evaluate the difference between references and results
    assert string_ref == string_result

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)
