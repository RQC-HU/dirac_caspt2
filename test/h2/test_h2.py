import subprocess
import os
import pytest
import glob


def test_h2():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path)  # Change directory to dirac_caspt2/test/h2
    print(test_path)  # Debug output
    ref_filename = "reference.H2.out"  # Reference
    ref_file_path = os.path.normpath(os.path.join(test_path, "/", ref_filename))
    output_filename = "H2.caspt2.out"  # Output (This file is compared with Reference)
    output_file_path = os.path.normpath(os.path.join(test_path, "/", output_filename))
    r4dcasci = os.path.normpath(os.path.join(test_path, "../../bin/r4dcascicoexe"))
    r4dcaspt2 = os.path.normpath(os.path.join(test_path, "../../bin/r4dcaspt2ocoexe"))
    # Run calculation
    p = subprocess.run(
        " ".join(
            [r4dcasci, "&>", output_file_path, "&&", r4dcaspt2, "&>>", output_file_path]
        ),
        shell=True,
    )
    print(
        "CASCI/CASPT2 status", p.returncode
    )  # Check status (If p.returncode != 0, calculation failed.)
    delete_files = [
        "[A-H]*int*",
        "MDCINTNEW",
        "NEWCICOEFF",
        "CIMAT*",
        "e0after",
        "EPS",
        "TRANSFOCK",
    ]

    # Delete scratch files
    for d in delete_files:
        files = glob.glob("".join([test_path, "/", d]))
        for f in files:
            os.remove(f)
    # Check output
    with open(output_file_path, encoding="utf-8", mode="r") as f:
        print(f.read())
    # Grep the reference output file
    with open(ref_file_path, encoding="utf-8", mode="r") as f:
        grep_str = [
            s.strip() for s in f.readlines() if "Total energy is" in s
        ]  # (e.g. ['Total energy is             -1.117672932144052 a.u.'])
        ref_energy = float(grep_str[-1].split()[-2])  # (e.g. -1.117672932144052)
    # Grep the test output file
    with open(output_file_path, encoding="utf-8", mode="r") as f:
        grep_str2 = [s.strip() for s in f.readlines() if "Total energy is" in s]
        output_energy = float(grep_str2[-1].split()[-2])

    # The previous output file is overwritten by the current output file.
    prev_output = "H2.caspt2.out.prev"  # Previous output (After test, the output file is moved to this)
    prev_file_path = os.path.normpath(os.path.join(test_path, "/", prev_output))
    subprocess.run(" ".join(["mv", output_filename, prev_file_path]), shell=True)
    # Check whether the output of test run
    # matches the reference to 7th decimal places.
    assert output_energy == pytest.approx(ref_energy, 1e-8)


if __name__ == "__main__":
    test_h2()
