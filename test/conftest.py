import pytest

slow_only_option = "--slowonly"
dev_option = "--dev"
runall_option = "--all"


def pytest_addoption(parser):
    parser.addoption(slow_only_option, action="store_true", default=False, help="run only very slow tests")
    parser.addoption(dev_option, action="store_true", default=False, help="run tests for development")
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
    config.addinivalue_line("markers", "slowonly: mark test as slow to run")
    config.addinivalue_line("markers", "dev: mark test as for development")


def pytest_collection_modifyitems(config, items):
    skip_slow = pytest.mark.skip(reason=f"need {runall_option} or {slow_only_option} option to run. REASON: slow test")
    skip_tests_because_dev = pytest.mark.skip(reason=f"need no option or {runall_option} option to run. REASON: --dev was activated")
    skip_fast_dev = pytest.mark.skip(reason=f"need no option or {dev_option} or {runall_option} option or  to run. REASON: --slowonly was activated")
    skip_fast_neutral = pytest.mark.skip(reason=f"need no option or {runall_option} option or  to run. REASON: --slowonly was activated")
    if config.getoption(runall_option):
        print("run all tests")
        return
    for item in items:
        # Check whether the test should be skipped or not.
        # The tests with item.add_marker("something") added will be skipped.

        # Tests marked by @pytest.mark.dev
        if item.get_closest_marker("dev"):
            # dev tests always run except when --slowonly was activated
            if config.getoption(slow_only_option):
                item.add_marker(skip_fast_dev)
            else:
                pass
        # Tests marked by @pytest.mark.slowonly
        elif item.get_closest_marker("slowonly"):
            # slow tests only run when --slowonly was activated
            if config.getoption(slow_only_option):
                pass
            else:
                item.add_marker(skip_slow)
        # Unmarked tests
        else:
            # Skip neutral tests if --dev or --slowonly were activated
            if config.getoption(dev_option):
                item.add_marker(skip_tests_because_dev)
            elif config.getoption(slow_only_option):
                item.add_marker(skip_fast_neutral)
