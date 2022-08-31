import pytest

slow_opion = "--slow"
very_slow_only_opion = "--veryslowonly"
runall_option = "--all"


def pytest_addoption(parser):
    parser.addoption(slow_opion, action="store_true", default=False, help="run normal tests and slow tests")
    parser.addoption(very_slow_only_opion, action="store_true", default=False, help="run only very slow tests")
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
    config.addinivalue_line("markers", "veryslowonly: run only very slow tests")


def pytest_collection_modifyitems(config, items):
    skip_slow = pytest.mark.skip(reason=f"need {runall_option} or {slow_opion} option to run. REASON: slow test")
    skip_fast_because_very_slow_only = pytest.mark.skip(reason=f"need no option or {runall_option} or {slow_opion} option to run. REASON: only very slow tests")
    skip_slow_because_very_slow_only = pytest.mark.skip(reason=f"need {runall_option} or {slow_opion} option to run. REASON: only very slow tests")
    skip_very_slow = pytest.mark.skip(reason=f"need {runall_option} or {very_slow_only_opion} option to run. REASON: very slow test")
    if config.getoption(runall_option):
        print("run all tests")
        return
    for item in items:
        if item.get_closest_marker("veryslowonly"):
            if not config.getoption(very_slow_only_opion):
                item.add_marker(skip_very_slow)
        elif item.get_closest_marker("slow"):
            if config.getoption(slow_opion):
                pass
            elif config.getoption(very_slow_only_opion):
                item.add_marker(skip_slow_because_very_slow_only)
            else:  # no args related to markers
                item.add_marker(skip_slow)
        else:  # Neutral tests
            if config.getoption(very_slow_only_opion):
                item.add_marker(skip_fast_because_very_slow_only)
