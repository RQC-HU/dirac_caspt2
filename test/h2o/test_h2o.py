import shutil
import subprocess
import os
import sys
import pytest
import glob


# Delete delete_files in the test_path
def delete_scratch_files(delete_files: "list[str]", test_path: str) -> None:
    for d in delete_files:
        files = glob.glob(os.path.abspath(os.path.join(test_path, d)))
        for f in files:
            os.remove(f)


def test_h2o(the_number_of_process: int) -> None:

    # Set file names
    ref_filename = "reference.H2O.out"  # Reference
    output_filename = "H2O.caspt2.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.H2O.caspt2.out"  # latest passed output (After test, the output file is moved to this)

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    binary_dir = os.path.abspath(
        os.path.join(test_path, "../../bin")
    )  # Set the Built binary directory
    r4dcasci: str = os.path.abspath(
        os.path.join(binary_dir, "r4dcascicoexe")
    )  # CASCI binary
    r4dcaspt2: str = os.path.abspath(
        os.path.join(binary_dir, "r4dcaspt2ocoexe")
    )  # CASPT2 binary

    # Set delete file list
    delete_files: list[str] = [
        "[A-H]*int*",  # 2-integrals per subspace
        "MDCINTNEW*",  # 2-integrals per MPI process
        "NEWCICOEFF",  # Coefficients of CI
        "CIMAT*",  # CI matrix
        "e0after",  # Energy after CASCI
        "EPS",  # epsilon
        "TRANSFOCK",  # Transformation matrix
        "*mat*",  # Matrix for DMRG
        "fort.*",  # Fortran files
    ]

    # Delete files because of previous test may illegally failed and created files that are not expected
    delete_scratch_files(delete_files, test_path)

    # Check binary files are exist
    if os.path.exists(r4dcasci) is False:
        error_message = (
            f"ERROR: {r4dcasci} is not exist.\nPlease build {r4dcasci} first."
        )
        print(error_message, file=sys.stderr)
        # Exit with error message
        sys.exit(error_message)
    if os.path.exists(r4dcaspt2) is False:
        error_message = (
            f"ERROR: {r4dcaspt2} is not exist.\nPlease build {r4dcaspt2} first."
        )
        print(error_message, file=sys.stderr)
        # Exit with error message
        sys.exit(error_message)

    # Set test command
    test_command = ""
    if the_number_of_process > 1:  # If the number of process is greater than 1, use MPI
        test_command = f"mpirun -np {the_number_of_process} {r4dcasci} && mpirun -np {the_number_of_process} {r4dcaspt2}"
    else:  # If the number of process is 1, use serial
        test_command = f"{r4dcasci} && {r4dcaspt2}"
    # Run calculation
    with open(output_file_path, "w") as file_output:
        p = subprocess.run(
            test_command,
            shell=True,
            encoding="utf-8",
            stdout=file_output,  # Redirect output to file_output
        )
    status = "CASCI/CASPT2 status " + str(p.returncode)
    # If the return code is not 0, print error message, probably calculation failed
    if p.returncode != 0:
        print(status, file=sys.stderr)

    # Delete scratch files
    delete_scratch_files(delete_files, test_path)

    # Check output
    with open(ref_file_path, encoding="utf-8", mode="r") as file_ref:
        try:  # Try to get the reference data
            # (e.g. ['Total energy is             -1.117672932144052 a.u.'])
            grep_str_ref: list[str] = [
                s.strip() for s in file_ref.readlines() if "Total energy is" in s
            ]
            ref_energy = float(
                grep_str_ref[-1].split()[-2]
            )  # (e.g. -1.117672932144052)
        except Exception as error:  # Failed to get the reference data
            error_message = f"{error}\nERROR: Failed to get the CASPT2 energy from the reference file {ref_file_path}."
            # Exit with error message
            sys.exit(error_message)

    # Grep the test output file
    with open(output_file_path, encoding="utf-8", mode="r") as file_output:
        try:  # Try to get the test data
            grep_str_output: list[str] = [
                s.strip() for s in file_output.readlines() if "Total energy is" in s
            ]
            output_energy = float(grep_str_output[-1].split()[-2])
        except Exception as error:  # Failed to get the test data
            error_message = f"{error}\nERROR: Failed to get the CASPT2 energy from the test file {output_file_path}."
            # Exit with error message
            sys.exit(error_message)

    # Check whether the output of test run
    # matches the reference to 7th decimal places.
    assert output_energy == pytest.approx(ref_energy, abs=1e-8)

    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)
