import subprocess
import os

def runtest():

    # Initialization
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Run calculation
    p = subprocess.run(
        " ".join(["pytest --tb=short -vv"]),
        shell=True,
        encoding="utf-8",
    )


if __name__ == "__main__":
    runtest()
