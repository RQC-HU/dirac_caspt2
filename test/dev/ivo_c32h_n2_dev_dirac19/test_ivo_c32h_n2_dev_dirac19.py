import os
import shutil

import pytest
from module_testing import (
    run_test_dcaspt2,
)


@pytest.mark.dev
def test_ivo_c32h_n2_dev_dirac19(env_setup_ivo) -> None:
    (test_path, DFPCMONEW_path, ref_DFPCMONEW_path, latest_passed_DFPCMONEW_path, output_path, latest_passed_output_path, test_command) = env_setup_ivo
    os.chdir(test_path)  # Change directory to the path of this file
    run_test_dcaspt2(test_command)

    # DFPCMONEW format
    # INFO (DIRAC version >= 21)
    # N2                                    Thu Apr   6 11:00:00 2023
    #  2 15 15 54 15 15 54
    #  -0.1075307794799569E+03
    # COEFS (DIRAC version >= 21)
    #     0.0783846162631894   -0.0932522358717089    0.2444662687107759   -0.2100050908725506    0.0207980763363816   -0.0061525045832165
    #    -0.0001106259309856   -0.0000860939270339    0.0001830653248163    0.0000000555844586   -0.0000000349239622    0.0000000146384444
    #    -0.0000000555844586    0.0000000349239622   -0.0000000146384444   -0.4656826540533246    0.2259566676583494   -0.3303017764820634
    # ...
    # EVALS (DIRAC version >= 21)
    #    -0.378466132952E+05   -0.376317160433E+05   -0.375871872904E+05   -0.375847489191E+05   -0.375740444395E+05   -0.375738772125E+05
    #    -0.375678183830E+05   -0.375638423228E+05   -0.375612685737E+05   -0.375611284085E+05   -0.375605808599E+05   -0.375591793895E+05
    #    -0.375586003412E+05   -0.375585970263E+05   -0.375581434458E+05   -0.155806087301E+02   -0.122487085753E+01   -0.527584459237E+00
    # ...
    # SUPERSYM
    #  1 1 1 1 -3 1 1 1 -3 1 1 1 -3 1 1 1 1 1 1 -3 1 1 -3 1 1 1 -3 1 1 1 1 1 1 1 1 -3 1 1 1 -3 1 1 1 -3 1 1 1 1 -3 1 1 -3 1 1 1 -3 1 1 1 1
    #  1 2 3 2

    # Open DFPCMONEW and reference.DFPCMONEW and compare the values (if the values are float, compare the values to 10th decimal places)
    with open(DFPCMONEW_path, "r") as DFPCMONEW_file, open(ref_DFPCMONEW_path, "r") as ref_file:
        for DFPCMONEW_line, ref_line in zip(DFPCMONEW_file, ref_file):
            # if the first value cannot be converted to float, compare the values as strings
            DFPCMONEW_values = DFPCMONEW_line.split()
            ref_values = ref_line.split()
            try:
                DFPCMONEW_float_values = [float(value) for value in DFPCMONEW_values]
                ref_float_values = [float(value) for value in ref_values]
                for DFPCMONEW_float_value, ref_float_value in zip(DFPCMONEW_float_values, ref_float_values):
                    assert ref_float_value == pytest.approx(DFPCMONEW_float_value, abs=1e-13)
            except ValueError:
                assert DFPCMONEW_values == ref_values

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_path, latest_passed_output_path)
    shutil.copy(DFPCMONEW_path, latest_passed_DFPCMONEW_path)
