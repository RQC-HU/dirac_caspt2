import subprocess
import os
import shutil
import sys
import pytest


def test_int_sort():

    # Set file names
    output_filename = "int.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.int.out"  # latest passed output (After test, the output file is moved to this)
    exe_filename = "test_sort_int_exe"  # Executable file

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    # Run tests
    p = subprocess.run(exe_file_path, shell=True)
    status = test_path + "status " + str(p.returncode)
    # If the return code is not 0, print error message, probably calculation failed
    if p.returncode != 0:
        print(status, file=sys.stderr)

    # Reference data
    reference_list:list[int] = [
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        169,
        170,
        171,
        172,
        173,
        174,
        175,
        1,
        3,
        5,
        156,
        189,
    ]
    reference_list.sort()  # 1,3,5,8,9,10,11,12,13,14,15,16,156,169,170,171,172,173,174,175,189

    # Get values from result
    with open(output_file_path) as file_result:
        try:  # Try to get the result data
            string_result = file_result.read()
            string_result = string_result.strip().split()
            result_real_list = list(map(float, string_result))
        except Exception as error:  # Failed to get the result data
            error_message = f"{error}\nERROR: Failed to get the data from the test file {output_file_path}."
            # Exit with error message
            sys.exit(error_message)

    # Evaluate the difference between references and results
    for out, ref in zip(result_real_list, reference_list):
        assert ref == out

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)


def test_int_sort_reverse():

    # Set file names
    output_filename = "int_reverse.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.int_reverse.out"  # latest passed output (After test, the output file is moved to this)
    exe_filename = "test_sort_int_reverse_exe"  # Executable file

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    # Run tests
    p = subprocess.run(exe_file_path, shell=True)
    status = test_path + "status " + str(p.returncode)
    # If the return code is not 0, print error message, probably calculation failed
    if p.returncode != 0:
        print(status, file=sys.stderr)

    # Reference data
    reference_list: list[int] = [
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        169,
        170,
        171,
        172,
        173,
        174,
        175,
        1,
        3,
        5,
        156,
        189,
    ]
    reference_list.sort(
        reverse=True
    )  # 189,175,174,173,172,171,170,169,156,16,15,14,13,12,11,10,9,8,5,3,1

    # Get values from result
    with open(output_file_path) as file_result:
        try:  # Try to get the result data
            string_result = file_result.read()
            string_result = string_result.strip().split()
            result_real_list = list(map(float, string_result))
        except Exception as error:  # Failed to get the result data
            error_message = f"{error}\nERROR: Failed to get the data from the test file {output_file_path}."
            # Exit with error message
            sys.exit(error_message)

    # Evaluate the difference between references and results
    for out, ref in zip(result_real_list, reference_list):
        assert ref == out

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)


def test_real_sort():

    # Set file names
    output_filename = "real.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.real.out"  # latest passed output (After test, the output file is moved to this)
    exe_filename = "test_sort_real_exe"  # Executable file

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    # Run tests
    p = subprocess.run(exe_file_path, shell=True)
    status = test_path + "status " + str(p.returncode)
    # If the return code is not 0, print error message, probably calculation failed
    if p.returncode != 0:
        print(status, file=sys.stderr)

    # Reference data
    reference_list: list[float] = [8.1, -9.2, 10000.58, -897, 123456789, 0.0000000010]
    reference_list.sort()  # -897, -9.2, 0.0000000010, 8.1, 10000.58, 123456789

    # Get values from result
    with open(output_file_path) as file_result:
        try:  # Try to get the result data
            string_result = file_result.read()
            string_result = string_result.strip().split()
            result_real_list = list(map(float, string_result))
        except Exception as error:  # Failed to get the result data
            error_message = f"{error}\nERROR: Failed to get the data from the test file {output_file_path}."
            # Exit with error message
            sys.exit(error_message)

    # Evaluate the difference between references and results
    for out, ref in zip(result_real_list, reference_list):
        assert ref == pytest.approx(out, 5e-7)

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)


def test_real_sort_reverse():

    # Set file names
    output_filename = (
        "real_reverse.out"  # Output (This file is compared with Reference)
    )
    latest_passed_output = "latest_passed.real_reverse.out"  # latest passed output (After test, the output file is moved to this)
    exe_filename = "test_sort_real_reverse_exe"  # Executable file

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    # Run tests
    p = subprocess.run(exe_file_path, shell=True)
    status = test_path + "status " + str(p.returncode)
    # If the return code is not 0, print error message, probably calculation failed
    if p.returncode != 0:
        print(status, file=sys.stderr)

    # Reference data
    reference_list: list[float] = [8.1, -9.2, 10000.58, -897, 123456789, 0.0000000010]
    reference_list.sort(
        reverse=True
    )  # 123456789, 10000.58, 8.1, 0.0000000010, -9.2, -897

    # Get values from result
    with open(output_file_path) as file_result:
        try:  # Try to get the result data
            string_result = file_result.read()
            string_result = string_result.strip().split()
            result_real_list: list[float] = list(map(float, string_result))
        except Exception as error:  # Failed to get the result data
            error_message = f"{error}\nERROR: Failed to get the data from the test file {output_file_path}."
            # Exit with error message
            sys.exit(error_message)

    # Evaluate the difference between references and results
    for out, ref in zip(result_real_list, reference_list):
        assert ref == pytest.approx(out, 5e-7)

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)
