import os
from module_testing import (
    create_test_command,
    run_test,
    get_split_string_list_from_output_file,
    is_binary_file_exist,
)
import pytest


@pytest.mark.dev
def test_ras3_bitcheck(env_setup_unittest):
    exe_file_path = env_setup_unittest("ras3_bitcheck_exe")
    # Current path
    test_path = os.path.dirname(os.path.abspath(__file__))

    # Change directory to the current path
    os.chdir(test_path)

    # input/output/executable file names
    ref_output_file = "expected"
    result_filename = "result.out"

    # Absolute path to input/output/executable files
    ref_output_file_path = os.path.abspath(os.path.join(test_path, ref_output_file))
    result_file_path = os.path.abspath(os.path.join(test_path, result_filename))

    is_binary_file_exist(exe_file_path)
    test_command = create_test_command(mpi_num_process=1, binaries=[exe_file_path])

    run_test(test_command)

    string_ref = get_split_string_list_from_output_file(ref_output_file_path)
    string_result = get_split_string_list_from_output_file(result_file_path)

    # Evaluate the difference between references and results
    assert string_ref == string_result


if __name__ == "__main__":
    test_ras3_bitcheck()
