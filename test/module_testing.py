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
        error_message = f"ERROR: {binary_file} is not exist.\nPlease build {binary_file} first."
        raise Exception(error_message)


def create_test_command_dcaspt2(dcaspt2: str, mpi_num_process: int, omp_num_threads: "int|None", input_file: str, output_file: str, test_path: str, save: bool, is_ivo: bool = False) -> str:
    options = ""
    if is_ivo:
        options = " --ivo --get \"DFPCMONEW\""
    if save:
        scratch_path = os.path.join(test_path, "scratch")
        options += f" --save --scratch {scratch_path}"
    if mpi_num_process > 1:
        options += f" --mpi {mpi_num_process}"
    if omp_num_threads is None:
        omp = str(os.environ.get("OMP_NUM_THREADS", 1))
    else:
        omp = str(omp_num_threads) if omp_num_threads > 1 else "1"
    options += f" --omp {omp}"

    test_command = f"{dcaspt2} -i {input_file} -o {output_file} {options}"
    return test_command


def create_test_command(mpi_num_process: int, binaries: "list[str]") -> str:
    test_command = ""
    if mpi_num_process > 1:  # If the number of process is greater than 1, use MPI
        for idx, binary in enumerate(binaries):
            if idx == 0:
                test_command = f"mpirun -np {mpi_num_process} {binary}"
            else:
                test_command = f"{test_command} && mpirun -np {mpi_num_process} {binary}"
    else:  # If the number of process is 1, use serial
        for idx, binary in enumerate(binaries):
            if idx == 0:
                test_command = f"{binary}"
            else:
                test_command = f"{test_command} && {binary}"
    with open("execution_command.log", "w") as file_output:
        file_output.write(test_command)
    return test_command


def run_test_dcaspt2(test_command: str) -> None:
    process = subprocess.run(test_command, shell=True)
    process.check_returncode()
    return


def run_test(test_command: str, output_file_path: str = "stdout.out") -> None:
    with open(output_file_path, "w") as file_output:
        process = subprocess.run(
            test_command,
            shell=True,
            encoding="utf-8",
            stdout=file_output,
        )
    process.check_returncode()
    return


def check_test_returncode(process: "subprocess.CompletedProcess[str]") -> None:
    if process.returncode != 0:
        raise Exception("ERROR: Process failed. return code status : " + str(process.returncode))


def get_caspt2_energy_from_output_file(file_path: str) -> float:
    with open(file_path, encoding="utf-8", mode="r") as output_file:
        try:
            # (e.g. ['Total energy is             -1.117672932144052 a.u.'])
            grep_str: list[str] = [s.strip() for s in output_file.readlines() if "Total energy is" in s]
            # (e.g. -1.117672932144052)
            caspt2_energy = float(grep_str[-1].split()[-2])
            return caspt2_energy
        except Exception as error:  # Failed to get the reference data
            error_message = f"{error}\nERROR: Failed to get the CASPT2 energy from the reference file {file_path}."
            raise Exception(error_message)


def get_stripped_string_from_output_file(file_path: str) -> str:
    with open(file_path, encoding="utf-8", mode="r") as output_file:
        try:
            string = output_file.read()
            return string.strip()
        except Exception as error:  # Failed to get the reference data
            error_message = f"{error}\nERROR: Failed to get the data from the reference file {file_path}."
            raise Exception(error_message)


def get_split_string_list_from_output_file(file_path: str) -> "list[str]":
    with open(file_path, encoding="utf-8", mode="r") as output_file:
        try:
            string = output_file.read()
            return string.strip().split()
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
