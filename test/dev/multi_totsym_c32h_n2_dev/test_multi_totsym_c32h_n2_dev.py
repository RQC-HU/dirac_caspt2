import os
import shutil
import pytest
from module_testing import (
    get_multi_caspt2_energy_from_output_file,
    run_test_dcaspt2,
)

@pytest.mark.dev
def test_multi_totsym_c32h_n2_dev(env_setup_caspt2) -> None:
    (test_path, ref_output_path, output_path, latest_passed_path, test_command) = env_setup_caspt2

    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    run_test_dcaspt2(test_command)

    ref_energies = get_multi_caspt2_energy_from_output_file(ref_output_path)
    test_energies = get_multi_caspt2_energy_from_output_file(output_path)

    assert len(ref_energies) == len(test_energies), "The number of energies in the reference and test output files do not match."

    # Check whether the output of test run
    # matches the reference to 10th decimal places.
    for idx, (ref_energy, test_energy) in enumerate(zip(ref_energies, test_energies)):
        assert test_energy == pytest.approx(ref_energy, abs=1e-10), f"Energy mismatch at index {idx}: {test_energy} != {ref_energy}"

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_path, latest_passed_path)
