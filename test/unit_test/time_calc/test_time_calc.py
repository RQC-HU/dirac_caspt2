import os
import shutil
from module_testing import (
    create_test_command,
    delete_scratch_files,
    get_stripped_string_from_output_file,
    get_split_string_list_from_output_file,
    is_binary_file_exist,
    run_test,
)
import pytest


@pytest.mark.dev
def test_time_calc():
    from datetime import datetime

    def create_date_time_from_input(input_str: str) -> datetime:
        # (e.g. "2021 01 01 00 00 00 00\n")
        splitted = input_str.rstrip().split()
        year_to_sec, ms_str = splitted[:-1], splitted[-1]
        year_to_sec[0] = str(int(year_to_sec[0])).zfill(4)  # Year must be 4 digits, but it might be less than 4 digits (e.g. 21 instead of 0021)
        microsecond = int(ms_str) * 1000  # input is in milliseconds, so times 1000 to covert to microseconds
        # year_to_sec + microsecond
        res_datetime = datetime.strptime(" ".join(year_to_sec), "%Y %m %d %H %M %S") + timedelta(microseconds=microsecond)
        return res_datetime

    def create_date_time_from_output(output_str: str) -> datetime:
        # (e.g. "computational time = 19 day 22 h 49 min 56.997 sec\n")
        splitted = output_str.rstrip().split()
        day, hour, minute, second = int(splitted[-8]), int(splitted[-6]), int(splitted[-4]), float(splitted[-2])
        return timedelta(days=day, hours=hour, minutes=minute, seconds=second)

    # Set file names
    input_file = "input"
    output_filename = "stdout.out"  # Output (This file is compared with Reference)
    latest_passed_output = "latest_passed.result.out"  # latest passed output (After test, the output file is moved to this)
    exe_filename = "test_time_calc_exe"  # Executable file

    # Get this files path and change directory to this path
    test_path = os.path.dirname(os.path.abspath(__file__))  # The path of this file
    os.chdir(test_path)  # Change directory to the path of this file
    print(test_path, "test start")  # Debug output

    # Set file paths
    input_file_path = os.path.abspath(os.path.join(test_path, input_file))
    output_file_path = os.path.abspath(os.path.join(test_path, output_filename))
    latest_passed_path = os.path.abspath(os.path.join(test_path, latest_passed_output))
    exe_file_path = os.path.abspath(os.path.join(test_path, exe_filename))

    is_binary_file_exist(exe_file_path)
    delete_scratch_files([output_filename], test_path)
    test_command = create_test_command(mpi_num_process=1, binaries=[exe_file_path])

    run_test(test_command)

    from datetime import datetime, timedelta
    from typing import List, Tuple

    # Create reference durations
    ref_durations: List[timedelta] = []
    start_and_end_times: List[Tuple[datetime, datetime]] = []
    with open(input_file_path, "r") as f:
        num_time_durations = int(f.readline().split()[0])
        for _ in range(num_time_durations):
            # (e.g. "2021 01 01 00 00 00 00")
            start_time = create_date_time_from_input(f.readline())
            end_time = create_date_time_from_input(f.readline())
            start_and_end_times.append((start_time, end_time))
            ref_durations.append(end_time - start_time)


    # Get result durations
    result_durations: List[timedelta] = []
    with open(output_file_path, "r") as f:
        for line in f:
            if "computational time" in line:
                result_durations.append(create_date_time_from_output(line))

    for ref, res, start_and_end in zip(ref_durations, result_durations, start_and_end_times):
        assert ref == res, f"Referecne != Result, Reference: {ref}, Result: {res}, \
Start time: {start_and_end[0]}, End time: {start_and_end[1]}"
    # If it reaches this point, the result of assert is true.
    # The latest passed output file is overwritten by the current output file if assert is True.
    shutil.copy(output_file_path, latest_passed_path)
