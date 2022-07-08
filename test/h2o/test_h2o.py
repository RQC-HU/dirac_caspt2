import shutil
import subprocess
import os
import sys
import pytest
import glob


def test_h2o():

    # Initialization
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output
    # Set file names
    ref_filename = "reference.H2O.out"  # Reference
    output_filename = "H2O.caspt2.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.H2O.caspt2.out"  # latest passed output (After test, the output file is moved to this)
    # Set file paths
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    binary_dir = os.path.abspath(
        os.path.join(test_path, "../../bin")
    )  # Set the Built binary directory
    r4dcasci = os.path.abspath(
        os.path.join(binary_dir, "r4dcascicoexe")
    )  # CASCI binary
    r4dcaspt2 = os.path.abspath(
        os.path.join(binary_dir, "r4dcaspt2ocoexe")
    )  # CASPT2 binary

    # Check binary files are exist
    if os.path.exists(r4dcasci) is False:
        error_message = (
            f"ERROR: {r4dcasci} is not exist.\nPlease build {r4dcasci} first."
        )
        print(error_message)
        # Exit with error message
        sys.exit(error_message)
    if os.path.exists(r4dcaspt2) is False:
        error_message = (
            f"ERROR: {r4dcaspt2} is not exist.\nPlease build {r4dcaspt2} first."
        )
        print(error_message)
        # Exit with error message
        sys.exit(error_message)

    # Run calculation
    with open(output_file_path, "w") as f:
        test_command = " ".join(["mpiexec -np 4", r4dcasci, "&&", "mpiexec -np 4", r4dcaspt2])
        p = subprocess.run(
            test_command,
            shell=True,
            encoding="utf-8",
            stdout=f,
            stderr=f,
        )
    status = "CASCI/CASPT2 status " + str(p.returncode)
    print(status)  # Print status (If p.returncode != 0, probably calculation failed.)

    # Set delete file list
    delete_files = [
        "[A-H]*int*",  # 2-integrals per subspace
        "MDCINTNEW*",  # 2-integrals per MPI process
        "NEWCICOEFF",  # Coefficients of CI
        "CIMAT*",  # CI matrix
        "e0after",  # Energy after CASCI
        "EPS",  # epsilon
        "TRANSFOCK",  # Transformation matrix
        "*mat*",  # Matrix for DMRG
    ]

    # Delete scratch files
    for d in delete_files:
        files = glob.glob(os.path.abspath(os.path.join(test_path, d)))
        for f in files:
            os.remove(f)

    # Check output
    with open(ref_file_path, encoding="utf-8", mode="r") as f:
        try:  # Try to get the reference data
            # (e.g. ['Total energy is             -1.117672932144052 a.u.'])
            grep_str = [s.strip() for s in f.readlines() if "Total energy is" in s]
            ref_energy = float(grep_str[-1].split()[-2])  # (e.g. -1.117672932144052)
        except:  # Failed to get the reference data
            error_message = f"ERROR: Failed to get the CASPT2 energy from the reference file {ref_file_path}."
            # Exit with error message
            sys.exit(error_message)

    # Grep the test output file
    with open(output_file_path, encoding="utf-8", mode="r") as f:
        try:  # Try to get the test data
            grep_str2 = [s.strip() for s in f.readlines() if "Total energy is" in s]
            output_energy = float(grep_str2[-1].split()[-2])
        except:  # Failed to get the test data
            error_message = f"ERROR: Failed to get the CASPT2 energy from the test file {output_file_path}."
            # Exit with error message
            sys.exit(error_message)

    # Check whether the output of test run
    # matches the reference to 7th decimal places.
    assert output_energy == pytest.approx(ref_energy, 1e-8)
    # The latest passed output file is overwritten by the current output file if assert equals True.
    shutil.copy(output_file_path, latest_passed_path)


if __name__ == "__main__":
    test_h2o()
