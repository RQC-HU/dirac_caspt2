import subprocess
import shutil
import os
import pytest


def test_int_sort():

    # Current path
    test_path = os.path.dirname(os.path.abspath(__file__))

    # Change directory to the current path
    os.chdir(test_path)

    # input/output/executable file names
    result_filename = "result_int"
    move_filename = "result_int.prev"
    exe_filename = "test_sort_int_exe"

    # Absolute path to input/output/executable files
    result_file_path = os.path.abspath(os.path.join(test_path, result_filename))
    move_file_path = os.path.abspath(os.path.join(test_path, move_filename))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    # Run tests
    subprocess.run(exe_file_path, shell=True)

    # Reference data
    ref = [
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
    ref.sort()  # 1,3,5,8,9,10,11,12,13,14,15,16,156,169,170,171,172,173,174,175,189

    # Get values from result
    with open(result_file_path) as file_result:
        string_result = file_result.read()
        string_result = string_result.split()
    result_int_list = list(map(int, string_result))

    # Move result files to move_file_path
    shutil.move(result_file_path, move_file_path)

    # Evaluate the difference between references and results
    assert ref == result_int_list


def test_int_sort_reverse():

    # Current path
    test_path = os.path.dirname(os.path.abspath(__file__))

    # Change directory to the current path
    os.chdir(test_path)

    # input/output/executable file names
    result_filename = "result_int_reverse"
    move_filename = "result_int_reverse.prev"
    exe_filename = "test_sort_int_reverse_exe"

    # Absolute path to input/output/executable files
    result_file_path = os.path.abspath(os.path.join(test_path, result_filename))
    move_file_path = os.path.abspath(os.path.join(test_path, move_filename))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    # Run tests
    subprocess.run(exe_file_path, shell=True)

    # Reference data
    ref = [
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
    ref.sort(
        reverse=True
    )  # 189,175,174,173,172,171,170,169,156,16,15,14,13,12,11,10,9,8,5,3,1

    # Get values from result
    with open(result_file_path) as file_result:
        string_result = file_result.read()
        string_result = string_result.split()
    result_int_list = list(map(int, string_result))

    # Move result files to move_file_path
    shutil.move(result_file_path, move_file_path)

    # Evaluate the difference between references and results
    assert ref == result_int_list


def test_real_sort():
    # Current path
    test_path = os.path.dirname(os.path.abspath(__file__))

    # Change directory to the current path
    os.chdir(test_path)

    # input/output/executable file names
    result_filename = "result_real"
    move_filename = "result_real.prev"
    exe_filename = "test_sort_real_exe"

    # Absolute path to input/output/executable files
    result_file_path = os.path.abspath(os.path.join(test_path, result_filename))
    move_file_path = os.path.abspath(os.path.join(test_path, move_filename))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    # Run tests
    subprocess.run(exe_file_path, shell=True)

    # Reference data
    ref = [8.1, -9.2, 10000.58, -897, 123456789, 0.0000000010]
    ref.sort()  # -897, -9.2, 0.0000000010, 8.1, 10000.58, 123456789

    # Get values from result
    with open(result_file_path) as file_result:
        string_result = file_result.read()
        string_result = string_result.split()
    result_real_list = list(map(float, string_result))

    # Move result files to move_file_path
    shutil.move(result_file_path, move_file_path)

    # Evaluate the difference between references and results
    for out, ref in zip(result_real_list, ref):
        assert ref == pytest.approx(out, 5e-7)


def test_real_sort_reverse():
    # Current path
    test_path = os.path.dirname(os.path.abspath(__file__))

    # Change directory to the current path
    os.chdir(test_path)

    # input/output/executable file names
    result_filename = "result_real_reverse"
    move_filename = "result_real_reverse.prev"
    exe_filename = "test_sort_real_reverse_exe"

    # Absolute path to input/output/executable files
    result_file_path = os.path.abspath(os.path.join(test_path, result_filename))
    move_file_path = os.path.abspath(os.path.join(test_path, move_filename))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    # Run tests
    subprocess.run(exe_file_path, shell=True)

     # Reference data
    ref = [8.1, -9.2, 10000.58, -897, 123456789, 0.0000000010]
    ref.sort(reverse=True)  # 123456789, 10000.58, 8.1, 0.0000000010, -9.2, -897

    # Get values from result
    with open(result_file_path) as file_result:
        string_result = file_result.read()
        string_result = string_result.split()
    result_real_list = list(map(float, string_result))

    # Move result files to move_file_path
    shutil.move(result_file_path, move_file_path)

    # Evaluate the difference between references and results
    for out, ref in zip(result_real_list, ref):
        assert ref == pytest.approx(out, 5e-7)


if __name__ == "__main__":
    test_int_sort()
    test_int_sort_reverse()
    test_real_sort()
    test_real_sort_reverse()
