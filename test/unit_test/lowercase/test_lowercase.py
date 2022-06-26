import glob
import subprocess
import os
import shutil


def test_lowercase():

    # Current path
    test_path = os.path.dirname(os.path.abspath(__file__))

    # Change directory to the current path
    os.chdir(test_path)

    # input/output/executable file names
    ref_filename = "expected"
    result_filename = "result"
    move_filename = "result.prev"
    exe_filename = "test_lowercase_exe"

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
        string_ref = string_ref.strip()

    # Get values from result
    with open(result_file_path) as file_result:
        string_result = file_result.read()
        string_result = string_result.strip()

    # Move result files to move_file_path
    shutil.move(result_file_path, move_file_path)

    # Evaluate the difference between references and results
    assert string_ref == string_result


if __name__ == "__main__":
    test_lowercase()
