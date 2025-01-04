import os
import subprocess

import pytest
from module_testing import (
    get_split_string_list_from_output_file,
)


@pytest.mark.dev
def test_pass_subspace_a_and_b(env_setup_gen_restart_file) -> None:
    (gen_restart_path, test_path, input_path, expected_path) = env_setup_gen_restart_file
    os.chdir(test_path)
    command = f"{gen_restart_path} {input_path}"
    subprocess.run(command.split(), check=True)

    string_ref = get_split_string_list_from_output_file(expected_path)
    string_result = get_split_string_list_from_output_file("caspt2_restart_2_1")

    for ref, result in zip(string_ref, string_result):
        assert ref == result


@pytest.mark.dev
def test_pass_full_subspaces(env_setup_gen_restart_file) -> None:
    (gen_restart_path, test_path, input_path, expected_path) = env_setup_gen_restart_file
    os.chdir(test_path)
    command = f"{gen_restart_path} {input_path}"
    subprocess.run(command.split(), check=True)

    string_ref = get_split_string_list_from_output_file(expected_path)
    string_result = get_split_string_list_from_output_file("caspt2_restart_33_1")

    for ref, result in zip(string_ref, string_result):
        assert ref == result


@pytest.mark.dev
def test_pass_in_the_middle_subspace_broken(env_setup_gen_restart_file) -> None:
    (gen_restart_path, test_path, input_path, expected_path) = env_setup_gen_restart_file
    os.chdir(test_path)
    command = f"{gen_restart_path} {input_path}"
    subprocess.run(command.split(), check=True)

    string_ref = get_split_string_list_from_output_file(expected_path)
    string_result = get_split_string_list_from_output_file("caspt2_restart_3_1")

    for ref, result in zip(string_ref, string_result):
        assert ref == result


@pytest.mark.dev
def test_pass_multi_ciroots(env_setup_gen_restart_file) -> None:
    (gen_restart_path, test_path, input_path, expected_path) = env_setup_gen_restart_file
    os.chdir(test_path)
    print(f"input_path: {input_path}")
    command = f"{gen_restart_path} {input_path}"
    subprocess.run(command.split(), check=True)

    for trailing_file_str in ["_2_1", "_2_2"]:
        string_ref = get_split_string_list_from_output_file(f"{expected_path}{trailing_file_str}")
        string_result = get_split_string_list_from_output_file(f"caspt2_restart{trailing_file_str}")

        for ref, result in zip(string_ref, string_result):
            assert ref == result
    # string_ref = get_split_string_list_from_output_file(expected_path)
    # string_result = get_split_string_list_from_output_file("caspt2_restart_4_1")

    # for ref, result in zip(string_ref, string_result):
    #     assert ref == result

@pytest.mark.dev
def test_no_output(env_setup_gen_restart_file) -> None:
    (gen_restart_path, test_path, input_path, expected_path) = env_setup_gen_restart_file
    os.chdir(test_path)
    command = f"{gen_restart_path} {input_path}"
    p = subprocess.run(command.split(), capture_output=True)
    assert p.returncode == 0 and p.stdout.decode("utf-8") == "" and p.stderr.decode("utf-8") == ""


@pytest.mark.dev
def test_fail_duplicated_energies(env_setup_gen_restart_file) -> None:
    (gen_restart_path, test_path, input_path, expected_path) = env_setup_gen_restart_file
    os.chdir(test_path)
    command = f"{gen_restart_path} {input_path}"
    p = subprocess.run(command.split(), capture_output=True)
    assert p.returncode != 0 and "ValueError" in p.stderr.decode("utf-8")
