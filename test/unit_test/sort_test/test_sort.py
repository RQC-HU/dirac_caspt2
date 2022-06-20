import subprocess
import shutil
import os
import pytest


def test_int_sort():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path)
    out_file_name = "intout"
    subprocess.run(
        os.path.abspath(os.path.join(test_path, "sort_test_int")), shell=True
    )

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
    out_file_path = os.path.abspath(os.path.join(test_path, out_file_name))
    with open(out_file_path) as f2:
        s2 = f2.read()
        s2 = s2.split()
    s2_int_list = list(map(int, s2))
    shutil.move(out_file_path, os.path.abspath(os.path.join(test_path, "intout_prev")))
    assert ref == s2_int_list


def test_int_sort_reverse():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path)
    out_file_name = "intoutreverse"
    subprocess.run(
        os.path.abspath(os.path.join(test_path, "sort_test_int_reverse")), shell=True
    )

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
    out_file_path = os.path.abspath(os.path.join(test_path, out_file_name))
    with open(out_file_path) as f2:
        s2 = f2.read()
        s2 = s2.split()
    s2_int_list = list(map(int, s2))
    shutil.move(
        out_file_path, os.path.abspath(os.path.join(test_path, "intoutreverse_prev"))
    )
    assert ref == s2_int_list


def test_real_sort():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path)
    out_file_name = "realout"
    subprocess.run(
        os.path.abspath(os.path.join(test_path, "sort_test_real")), shell=True
    )

    ref = [8.1, -9.2, 10000.58, -897, 123456789, 0.0000000010]
    ref.sort()  # -897, -9.2, 0.0000000010, 8.1, 10000.58, 123456789
    out_file_path = os.path.abspath(os.path.join(test_path, out_file_name))
    with open(out_file_path) as f2:
        s2 = f2.read()
        s2 = s2.split()
    s2_int_list = list(map(float, s2))
    shutil.move(out_file_path, os.path.abspath(os.path.join(test_path, "realout_prev")))
    for out, ref in zip(s2_int_list, ref):
        assert ref == pytest.approx(out, 5e-7)


def test_real_sort_reverse():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path)
    out_file_name = "realoutreverse"
    subprocess.run(
        os.path.abspath(os.path.join(test_path, "sort_test_real_reverse")), shell=True
    )

    ref = [8.1, -9.2, 10000.58, -897, 123456789, 0.0000000010]
    ref.sort(reverse=True)  # 123456789, 10000.58, 8.1, 0.0000000010, -9.2, -897
    out_file_path = os.path.abspath(os.path.join(test_path, out_file_name))
    with open(out_file_path) as f2:
        s2 = f2.read()
        s2 = s2.split()
    s2_int_list = list(map(float, s2))
    shutil.move(
        out_file_path, os.path.abspath(os.path.join(test_path, "realoutreverse_prev"))
    )
    for out, ref in zip(s2_int_list, ref):
        assert ref == pytest.approx(out, 5e-7)


if __name__ == "__main__":
    test_int_sort()
    test_real_sort()
