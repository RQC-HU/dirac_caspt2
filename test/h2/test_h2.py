import subprocess
import os
import pytest


def test_h2():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path) # Change directory to
    output_filename = "H2.caspt2.out"
    p = subprocess.run("sh " + test_path + "/h2.sh", shell=True)  # Run calculation
    print(p.returncode)
    ref = subprocess.run(
        "cat "
        + test_path
        + "/reference.H2.out"
        + " |  awk ' /Total/{print}' | tr -s '\n'",
        shell=True,
        encoding="utf-8",
        stdout=subprocess.PIPE,
    ) # Get CASPT2 energy and e2sum from the output
    output = subprocess.run(
        "cat "
        + test_path
        + "/"
        + output_filename
        + " | awk '/Total/{print}' | tr -s '\n'",
        shell=True,
        encoding="utf-8",
        stdout=subprocess.PIPE,
    ) # Get CASPT2 energy and e2sum from the output
    # Convert string to float
    ref_energy = float(ref.stdout.split()[-2])
    output_energy = float(output.stdout.split()[-2])
    # Check whether the output of test run matches the reference to 7th decimal places.
    assert output_energy == pytest.approx(ref_energy, 1e-8)


if __name__ == "__main__":
    test_h2()
