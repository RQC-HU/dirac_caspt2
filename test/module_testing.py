import glob
import os
import subprocess


def delete_scratch_files(delete_files: "list[str]", test_path: str) -> None:
    for d in delete_files:
        files = glob.glob(os.path.abspath(os.path.join(test_path, d)))
        for f in files:
            os.remove(f)


def is_binary_file_exist(binary_file: str) -> None:
    if not os.path.exists(binary_file):
        error_message = (
            f"ERROR: {binary_file} is not exist.\nPlease build {binary_file} first."
        )
        raise Exception(error_message)


def create_test_command(the_number_of_process: int, binaries: "list[str]") -> str:
    if the_number_of_process > 1:  # If the number of process is greater than 1, use MPI
        for idx, binary in enumerate(binaries):
            if idx == 0:
                test_command = f"mpirun -np {the_number_of_process} {binary}"
            else:
                test_command += f" && mpirun -np {the_number_of_process} {binary}"
    else:  # If the number of process is 1, use serial
        for idx, binary in enumerate(binaries):
            if idx == 0:
                test_command = f"{binary}"
            else:
                test_command += f" && {binary}"
    return test_command


def run_test(
    test_command: str, output_file_path: str
) -> "subprocess.CompletedProcess[str]":
    with open(output_file_path, "w") as file_output:
        process = subprocess.run(
            test_command,
            shell=True,
            encoding="utf-8",
            stdout=file_output,  # Redirect output to file_output
        )
    return process


def check_test_returncode(process: "subprocess.CompletedProcess[str]") -> None:
    if process.returncode != 0:
        raise Exception(
            "ERROR: Process failed. return code status : " + str(process.returncode)
        )


def get_caspt2_energy_from_output_file(file_path: str) -> float:
    with open(file_path, encoding="utf-8", mode="r") as output_file:
        try:
            # (e.g. ['Total energy is             -1.117672932144052 a.u.'])
            grep_str: list[str] = [
                s.strip() for s in output_file.readlines() if "Total energy is" in s
            ]
            caspt2_energy = float(grep_str[-1].split()[-2])  # (e.g. -1.117672932144052)
            return caspt2_energy
        except Exception as error:  # Failed to get the reference data
            error_message = f"{error}\nERROR: Failed to get the CASPT2 energy from the reference file {file_path}."
            raise Exception(error_message)


def get_stripped_string_from_output_file(file_path: str) -> str:
    with open(file_path, encoding="utf-8", mode="r") as output_file:
        try:
            string = output_file.read()
            string = string.strip()
            return string
        except Exception as error:  # Failed to get the reference data
            error_message = f"{error}\nERROR: Failed to get the data from the reference file {file_path}."
            raise Exception(error_message)


def get_split_string_list_from_output_file(file_path: str) -> "list[str]":
    with open(file_path, encoding="utf-8", mode="r") as output_file:
        try:
            string = output_file.read()
            string = string.strip().split()
            return string
        except Exception as error:  # Failed to get the reference data
            error_message = f"{error}\nERROR: Failed to get the data from the reference file {file_path}."
            raise Exception(error_message)


def convert_string_list_to_integer_list(string_list: "list[str]") -> "list[int]":
    try:
        integer_list = list(map(int, string_list))
        return integer_list
    except Exception as error:  # Failed to get the reference data
        error_message = f"{error}\nERROR: Failed to convert the string list to integer list. string_list : {string_list}"
        raise Exception(error_message)


def convert_string_list_to_float_list(string_list: "list[str]") -> "list[float]":
    try:
        float_list = list(map(float, string_list))
        return float_list
    except Exception as error:  # Failed to get the reference data
        error_message = f"{error}\nERROR: Failed to convert the string list to float list. string_list : {string_list}"
        raise Exception(error_message)
