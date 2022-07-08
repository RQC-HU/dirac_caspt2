import shutil
import subprocess
import os
import pytest
import glob


def test_lower_MPI_h2():
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path)  # Debug output
    ref_filename = "reference.H2.out"  # Reference
    ref_file_path = os.path.abspath(os.path.join(test_path, ref_filename))
    output_filename = "H2.caspt2.out"  # Output (This file is compared with Reference)
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    bindir = os.path.abspath(
        os.path.join(test_path, "../../bin")
    )  # Build binary directory
    r4dcasci = os.path.abspath(os.path.join(bindir, "r4dcascicoexe"))  # CASCI
    r4dcaspt2 = os.path.abspath(os.path.join(bindir, "r4dcaspt2ocoexe"))  # CASPT2
    # Run calculation
    with open(output_file_path, "w") as f:
        p = subprocess.run(
            " ".join([r4dcasci, "&&", r4dcaspt2]),
            shell=True,
            encoding="utf-8",
            stdout=f,
            stderr=f,
        )
    f.close()
    print(
        "CASCI/CASPT2 status", p.returncode
    )  # Debug output, Check status (If p.returncode != 0, calculation failed.)
    delete_files = [
        "[A-H]*int*",
        "MDCINTNEW*",
        "NEWCICOEFF",
        "CIMAT*",
        "e0after",
        "EPS",
        "TRANSFOCK",
    ]

    # Delete scratch files
    for d in delete_files:
        files = glob.glob(os.path.abspath(os.path.join(test_path, d)))
        print("files", files)  # Debug output
        for f in files:
            os.remove(f)

    # Check output
    with open(ref_file_path, encoding="utf-8", mode="r") as f:
        # (e.g. ['Total energy is             -1.117672932144052 a.u.'])
        grep_str = [s.strip() for s in f.readlines() if "Total energy is" in s]
        ref_energy = float(grep_str[-1].split()[-2])  # (e.g. -1.117672932144052)
    f.close()
    # Grep the test output file
    with open(output_file_path, encoding="utf-8", mode="r") as f:
        grep_str2 = [s.strip() for s in f.readlines() if "Total energy is" in s]
        output_energy = float(grep_str2[-1].split()[-2])
    f.close()

    # The previous output file is overwritten by the current output file.
    prev_output = "H2.caspt2.out.prev"  # Previous output (After test, the output file is moved to this)
    prev_file_path = os.path.abspath(os.path.join(test_path, prev_output))
    shutil.move(output_file_path, prev_file_path)

    # Check whether the output of test run
    # matches the reference to 7th decimal places.
    assert output_energy == pytest.approx(ref_energy, 1e-8)


if __name__ == "__main__":
    test_lower_MPI_h2()
