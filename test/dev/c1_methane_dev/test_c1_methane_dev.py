import os
import shutil

import pytest
from module_testing import (
    get_caspt2_energy_from_output_file,
    run_test_dcaspt2,
)


@pytest.mark.dev
def test_c1_methane_dev(env_setup_caspt2) -> None:
    (test_path, ref_output_path, output_path, latest_passed_path, test_command) = env_setup_caspt2

    # Get this files path and change directory to this path
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    with open("execution_command.txt", "w") as f:
        print(f"TEST COMMAND: {test_command}", file=f)
    run_test_dcaspt2(test_command)

    ref_energy = get_caspt2_energy_from_output_file(ref_output_path)
    test_energy = get_caspt2_energy_from_output_file(output_path)

    # Check whether the output of test run
    # matches the reference to 7th decimal places.
    assert test_energy == pytest.approx(ref_energy, abs=1e-9)

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_path, latest_passed_path)
