import glob
import subprocess
import os


def test_input_reader():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path)
    ref_filename = "ref.ras3"
    subprocess.run(
        os.path.abspath(os.path.join(test_path, "ras_input_reader_exe")), shell=True
    )

    with open(ref_filename) as f:
        s = f.read()
        s = s.split()
    s_int_list = list(map(int, s))
    s_int_list.sort()
    with open("file") as f2:
        s2 = f2.read()
        s2 = s2.split()

    s2_int_list = list(map(int, s2))
    delete_files = ["file"]
    for d in delete_files:
        files = glob.glob(os.path.abspath(os.path.join(test_path, d)))
        for f in files:
            os.remove(f)
    assert s_int_list == s2_int_list


if __name__ == "__main__":
    test_input_reader()
