import pytest

slow_opion = "--runslow"
runall_option = "--runall"


def pytest_addoption(parser):
    parser.addoption(slow_opion, action="store_true", default=False, help="run slow tests")
    parser.addoption(runall_option, action="store_true", default=False, help="run all tests")
    parser.addoption(
        "--parallel",
        type=int,
        default=1,
        help="run tests in parallel processes",
    )


@pytest.fixture
def the_number_of_process(request):
    return request.config.getoption("--parallel")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")
    config.addinivalue_line("markers", "runall: mark test as runall to run")


def pytest_collection_modifyitems(config, items):
    skip_slow = pytest.mark.skip(reason=f"need {slow_opion} or {runall_option} option to run. REASON: Slow test")
    skip_slowest = pytest.mark.skip(reason=f"need {runall_option} option to run. REASON: Slow test")
    for item in items:
        if item.get_closest_marker("runall") and not config.getoption(runall_option):
            item.add_marker(skip_slowest)
        elif item.get_closest_marker("slow") and not config.getoption(slow_opion) and not config.getoption(runall_option):
            item.add_marker(skip_slow)
