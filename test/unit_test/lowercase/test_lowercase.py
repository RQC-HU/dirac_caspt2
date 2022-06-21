import glob
import subprocess
import os
import shutil


def test_lowercase():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path)
    ref_filename = "ref"
    subprocess.run(
        os.path.abspath(os.path.join(test_path, "test_lowercase_exe")), shell=True
    )

    with open(ref_filename) as f:
        s = f.read()
        s = s.strip()
    with open("output") as f2:
        s2 = f2.read()
        s2 = s2.strip()
    delete_file = "output"
    mv_file = "output.prev"
    file = os.path.abspath(os.path.join(test_path, delete_file))
    mv_file_path = os.path.abspath(os.path.join(test_path, mv_file))
    shutil.move(file, mv_file_path)
    assert s == s2


if __name__ == "__main__":
    test_lowercase()
