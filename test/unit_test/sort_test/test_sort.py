import os
import pytest
import shutil
from module_testing import (
    check_test_returncode,
    convert_string_list_to_float_list,
    convert_string_list_to_integer_list,
    create_test_command,
    delete_scratch_files,
    is_binary_file_exist,
    run_test,
    get_split_string_list_from_output_file,
)


@pytest.mark.dev
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

    is_binary_file_exist(exe_file_path)
    delete_scratch_files([output_filename], test_path)
    test_command = create_test_command(the_number_of_process=1, binaries=[exe_file_path])

    process = run_test(test_command)
    check_test_returncode(process)

    # Reference data
    reference_list = [
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

    string_result = get_split_string_list_from_output_file(output_file_path)
    result_int_list = convert_string_list_to_integer_list(string_result)

    # Evaluate the difference between references and results
    for out, ref in zip(result_int_list, reference_list):
        assert ref == out

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)


@pytest.mark.dev
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

    is_binary_file_exist(exe_file_path)
    delete_scratch_files([output_filename], test_path)
    test_command = create_test_command(the_number_of_process=1, binaries=[exe_file_path])

    process = run_test(test_command, output_file_path)
    check_test_returncode(process)

    # Reference data
    reference_list = [
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
    reference_list.sort(reverse=True)  # 189,175,174,173,172,171,170,169,156,16,15,14,13,12,11,10,9,8,5,3,1

    string_result = get_split_string_list_from_output_file(output_file_path)
    result_int_list = convert_string_list_to_integer_list(string_result)

    # Evaluate the difference between references and results
    for out, ref in zip(result_int_list, reference_list):
        assert ref == out

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)


@pytest.mark.dev
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

    is_binary_file_exist(exe_file_path)
    delete_scratch_files([output_filename], test_path)
    test_command = create_test_command(the_number_of_process=1, binaries=[exe_file_path])

    process = run_test(test_command, output_file_path)
    check_test_returncode(process)

    # Reference data
    reference_list: list[float] = [8.1, -9.2, 10000.58, -897, 123456789, 0.0000000010]
    reference_list.sort()  # -897, -9.2, 0.0000000010, 8.1, 10000.58, 123456789

    string_result = get_split_string_list_from_output_file(output_file_path)
    result_real_list = convert_string_list_to_float_list(string_result)

    # Evaluate the difference between references and results
    for out, ref in zip(result_real_list, reference_list):
        assert ref == pytest.approx(out)

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)


@pytest.mark.dev
def test_real_sort_reverse():

    # Set file names
    output_filename = "real_reverse.out"  # Output (This file is compared with Reference)
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

    is_binary_file_exist(exe_file_path)
    delete_scratch_files([output_filename], test_path)
    test_command = create_test_command(the_number_of_process=1, binaries=[exe_file_path])

    process = run_test(test_command, output_file_path)
    check_test_returncode(process)

    # Reference data
    reference_list: list[float] = [8.1, -9.2, 10000.58, -897, 123456789, 0.0000000010]
    reference_list.sort(reverse=True)  # 123456789, 10000.58, 8.1, 0.0000000010, -9.2, -897

    string_result = get_split_string_list_from_output_file(output_file_path)
    result_real_list = convert_string_list_to_float_list(string_result)

    # Evaluate the difference between references and results
    for out, ref in zip(result_real_list, reference_list):
        assert ref == pytest.approx(out)

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)
