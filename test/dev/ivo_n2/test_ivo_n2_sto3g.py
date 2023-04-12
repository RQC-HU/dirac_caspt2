import os
import shutil
import pytest
from module_testing import (
    run_test_dcaspt2,
    create_test_command_dcaspt2,
)


def run_ivo_test_n2_sto3g(mpi_num_process: int, omp_num_threads: int, save: bool, dirac_version: int, path_this_file: str, test_path: str) -> None:

    # Set file names
    input_file = "active.ivo.inp"  # Input
    DFPCMONEW_file = "DFPCMONEW"  # Test (This file is compared with Reference)
    ref_DFPCMONEW_file = f"reference.DFPCMONEW.dirac{dirac_version}"  # Reference
    latest_passed_test = f"latest_passed.DFPCMONEW.dirac{dirac_version}"  # latest passed DFPCMONEW
    output_filename = f"c32h_n2_dev.dirac{dirac_version}.caspt2.out"  # Output
    latest_passed_output = f"latest_passed.c32h_n2_dev.dirac{dirac_version}.caspt2.out"  # latest passed output (After test, the output file is moved to this)

    # Set file paths
    input_file_path = os.path.abspath(os.path.join(test_path, input_file))
    DFPCMONEW_file_path = os.path.abspath(os.path.join(test_path, DFPCMONEW_file))
    ref_DFPCMONEW_file_path = os.path.abspath(os.path.join(path_this_file, ref_DFPCMONEW_file))
    latest_passed_DFPCMONEW_file_path = os.path.abspath(os.path.join(path_this_file, latest_passed_test))
    output_file_path = os.path.abspath(os.path.join(path_this_file, output_filename))
    latest_passed_output_file_path = os.path.abspath(os.path.join(path_this_file, latest_passed_output))
    binary_dir = os.path.abspath(os.path.join(path_this_file, "../../../bin"))  # Set the Built binary directory
    dcaspt2 = os.path.join(binary_dir, "dcaspt2")  # Set the dcaspt2 binary path

    is_ivo = True
    test_command = create_test_command_dcaspt2(dcaspt2, mpi_num_process, omp_num_threads, input_file_path, output_file_path, test_path, save, is_ivo)
    run_test_dcaspt2(test_command)

    # Copy DFPCMO to path_this_file/result.DFPCMONEW.dirac{dirac_version}
    DFPCMONEW_copy_path = os.path.abspath(os.path.join(path_this_file, f"result.DFPCMONEW.dirac{dirac_version}"))
    shutil.copyfile(DFPCMONEW_file_path, DFPCMONEW_copy_path)

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
    with open(DFPCMONEW_file_path, "r") as DFPCMONEW_file, open(ref_DFPCMONEW_file_path, "r") as ref_file:
        for DFPCMONEW_line, ref_line in zip(DFPCMONEW_file, ref_file):
            # if the first value cannot be converted to float, compare the values as strings
            DFPCMONEW_values = DFPCMONEW_line.split()
            ref_values = ref_line.split()
            try:
                DFPCMONEW_float_values = [float(value) for value in DFPCMONEW_values]
                ref_float_values = [float(value) for value in ref_values]
                for DFPCMONEW_float_value, ref_float_value in zip(DFPCMONEW_float_values, ref_float_values):
                    pytest.approx(DFPCMONEW_float_value, ref_float_value, abs=1e-10)
            except ValueError:
                assert DFPCMONEW_values == ref_values

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_output_file_path)
    shutil.copy(DFPCMONEW_file_path, latest_passed_DFPCMONEW_file_path)


@pytest.mark.dev
def test_ivo_n2_sto3g(mpi_num_process: int, omp_num_threads: int, save: bool) -> None:

    for dirac_version in [19, 22]:
        # Get this files path and change directory to this path
        path_this_file = os.path.dirname(os.path.abspath(__file__))  # The path of this file
        test_path = os.path.join(path_this_file, f"./input/{dirac_version}")  # The path of the test directory
        os.chdir(test_path)

        run_ivo_test_n2_sto3g(mpi_num_process, omp_num_threads, save, dirac_version, path_this_file, test_path)
