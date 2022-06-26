import glob
import shutil
import subprocess
import os


def test_ras_input_reader():

    # Current path
    test_path = os.path.dirname(os.path.abspath(__file__))

    # Change directory to the current path
    os.chdir(test_path)

    # input/output/executable file names
    ref_filename = "expected"
    result_filename = "result"
    move_filename = "result.prev"
    exe_filename = "test_ras_input_reader_exe"

    # Absolute path to input/output/executable files
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    result_file_path = os.path.abspath(os.path.join(test_path, result_filename))
    move_file_path = os.path.abspath(os.path.join(test_path, move_filename))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    # Run tests
    subprocess.run(exe_file_path, shell=True)

    # Get values from reference
    with open(ref_file_path) as file_ref:
        string_ref = file_ref.read()
        string_ref = string_ref.split()
    ref_int_list = list(map(int, string_ref))

    # Get values from result
    with open(result_file_path) as file_result:
        string_result = file_result.read()
        string_result = string_result.split()
    result_int_list = list(map(int, string_result))

    # Move result files to move_file_path
    shutil.move(result_file_path, move_file_path)

    # Evaluate the difference between references and results
    assert ref_int_list == result_int_list


if __name__ == "__main__":
    test_ras_input_reader()
