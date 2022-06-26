import glob
import shutil
import subprocess
import os


def test_ras3_bitcheck():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path)
    ref_filename = "expected"
    output_filename = "result"
    move_filename = "result.prev"
    exe_filename = "ras3_bitcheck_exe"
    subprocess.run(os.path.abspath(os.path.join(test_path, exe_filename)), shell=True)

    with open(ref_filename) as file_ref:
        string_ref = file_ref.read()
        string_ref = string_ref.split()
    with open(output_filename) as file_result:
        string_result = file_result.read()
        string_result = string_result.split()

    shutil.move(
        os.path.abspath(os.path.join(test_path, output_filename)),
        os.path.abspath(os.path.join(test_path, move_filename)),
    )
    assert string_ref == string_result


if __name__ == "__main__":
    test_ras3_bitcheck()
