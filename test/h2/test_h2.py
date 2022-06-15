import subprocess
import os
import pytest


def test_h2():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path)  # Change directory to dirac_caspt2/test/h2
    print(test_path)
    ref_filename = "reference.H2.out"   # Reference
    output_filename = "H2.caspt2.out"   # Output (This file is compared with Reference)
    prev_output = "H2.caspt2.out.prev"  # Previous output (After test, the output file is moved to this)
    # Run calculation
    p = subprocess.run("sh ./h2.sh", shell=True)
    print("CASCI/CASPT2 status", p.returncode)  # Check status (If p.returncode != 0, calculation failed.)
    ref = subprocess.run(
        "cat " + ref_filename + " |  awk ' /Total/{print}' | tr -s '\n'",
        shell=True,
        encoding="utf-8",
        stdout=subprocess.PIPE,
    )  # Get CASPT2 energy and e2sum from the output
    output = subprocess.run(
        "cat " + output_filename + " | awk '/Total/{print}' | tr -s '\n'",
        shell=True,
        encoding="utf-8",
        stdout=subprocess.PIPE,
    )  # Get CASPT2 energy and e2sum from the output
    # Check output for debugging
    print("Check output for debugging.")
    print("ref.stdout", ref.stdout)
    print("ouput.stdout", output.stdout)
    # Convert string to float
    ref_energy = float(ref.stdout.split()[-2])
    output_energy = float(output.stdout.split()[-2])
    # The previous output file is overwritten by the current output file.
    subprocess.run("mv " + output_filename + " " + prev_output, shell=True)
    # Check whether the output of test run
    # matches the reference to 7th decimal places.
    assert output_energy == pytest.approx(ref_energy, 1e-8)


if __name__ == "__main__":
    test_h2()
