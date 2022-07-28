import shutil
import os
from module_testing import (
    check_test_returncode,
    create_test_command,
    run_test,
    get_split_string_list_from_output_file,
)


def test_ras3_bitcheck():
    # Current path
    test_path = os.path.dirname(os.path.abspath(__file__))

    # Change directory to the current path
    os.chdir(test_path)

    # input/output/executable file names
    ref_filename = "expected"
    result_filename = "result"
    move_filename = "result.prev"
    exe_filename = "ras3_bitcheck_exe"

    # Absolute path to input/output/executable files
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    result_file_path = os.path.abspath(os.path.join(test_path, result_filename))
    move_file_path = os.path.abspath(os.path.join(test_path, move_filename))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    test_command = create_test_command(
        the_number_of_process=1, binaries=[exe_file_path]
    )

    process = run_test(test_command, result_file_path)
    check_test_returncode(process)

    string_ref = get_split_string_list_from_output_file(ref_file_path)
    string_result = get_split_string_list_from_output_file(result_file_path)

    # Move result files to move_file_path
    shutil.move(result_file_path, move_file_path)

    # Evaluate the difference between references and results
    assert string_ref == string_result


if __name__ == "__main__":
    test_ras3_bitcheck()
