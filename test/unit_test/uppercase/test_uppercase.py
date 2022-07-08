import subprocess
import os
import shutil
import sys
import pytest

def test_uppercase():

    # Set file names
    ref_filename = "expected"  # Reference
    output_filename = "result.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.result.out"  # latest passed output (After test, the output file is moved to this)
    exe_filename = "test_uppercase_exe" # Executable file

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))


    # Run tests
    p = subprocess.run(exe_file_path, shell=True)
    status = test_path + "status " + str(p.returncode)
    # If the return code is not 0, print error message, probably calculation failed
    if p.returncode != 0:
        print(status, file=sys.stderr)


    # Get values from reference
    with open(ref_file_path) as file_ref:
        try:
            string_ref = file_ref.read()
            string_ref = string_ref.strip()
        except Exception as error:  # Failed to get the reference data
            error_message = f"{error}\nERROR: Failed to get the data from the reference file {ref_file_path}."
            # Exit with error message
            sys.exit(error_message)

    # Get values from result
    with open(output_file_path) as file_result:
        try: # Try to get the result data
            string_result = file_result.read()
            string_result = string_result.strip()
        except Exception as error:  # Failed to get the result data
            error_message = f"{error}\nERROR: Failed to get the data from the test file {output_file_path}."
            # Exit with error message
            sys.exit(error_message)

    # Evaluate the difference between references and results
    assert string_ref == string_result

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)
