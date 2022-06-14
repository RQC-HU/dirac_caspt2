import subprocess
import os
import pytest


def test_h2():
    test_path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(test_path)  # Change directory to dirac_caspt2/test/h2
    r4dcascicoexe = os.path.join(test_path, "../../bin/r4dcascicoexe")
    r4dcaspt2ocoexe = os.path.join(test_path, "../../bin/r4dcaspt2ocoexe")
    print(test_path, r4dcascicoexe, r4dcaspt2ocoexe)
    ref_filename = "reference.H2.out"
    ref_filepath = os.path.join(test_path, ref_filename)
    output_filename = "H2.caspt2.out"
    output_filepath = os.path.join(test_path, output_filename)
    # Run calculation
    # p = subprocess.run("sh " + test_path + "/h2.sh", shell=True)
    # CASCI
    p = subprocess.run(r4dcascicoexe + " &> " + output_filename, shell=True)
    print("CASCI status", p.returncode)
    # CASPT2
    p = subprocess.run(r4dcaspt2ocoexe + " &>> " + output_filename, shell=True)
    print("CASPT2 status", p.returncode)
    # Remove scratch files
    p = subprocess.run(
        "rm [A-H]*int* CIMAT* e0after EPS MDCINTNEW NEWCICOEFF TRANSFOCK", shell=True
    )
    p = subprocess.run("ls", shell=True)
    ref = subprocess.run(
        "cat " + ref_filepath + " |  awk ' /Total/{print}' | tr -s '\n'",
        shell=True,
        encoding="utf-8",
        stdout=subprocess.PIPE,
    )  # Get CASPT2 energy and e2sum from the output
    print("ref output awk end")
    output = subprocess.run(
        "cat " + output_filepath + " | awk '/Total/{print}' | tr -s '\n'",
        shell=True,
        encoding="utf-8",
        stdout=subprocess.PIPE,
    )  # Get CASPT2 energy and e2sum from the output
    print("ref output awk end")
    # Check output for debugging
    print("ref.stdout", ref.stdout)
    print("ref.stdout", output.stdout)
    # Convert string to float
    ref_energy = float(ref.stdout.split()[-2])
    output_energy = float(output.stdout.split()[-2])
    # Check whether the output of test run
    # matches the reference to 7th decimal places.
    assert output_energy == pytest.approx(ref_energy, 1e-8)


if __name__ == "__main__":
    test_h2()
