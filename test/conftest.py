import pytest

slow_opion = "--runslow"
slowest_opion = "--runslowest"


def pytest_addoption(parser):
    parser.addoption(slow_opion, action="store_true", default=False, help="run slow tests")
    parser.addoption(slowest_opion, action="store_true", default=False, help="run slow tests")
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
    config.addinivalue_line("markers", "slowest: mark test as slowest to run")


def pytest_collection_modifyitems(config, items):
    skip_slow = pytest.mark.skip(reason=f"need {slow_opion} or {slowest_opion} option to run. REASON: Slow test")
    skip_slowest = pytest.mark.skip(reason=f"need {slowest_opion} option to run. REASON: Slow test")
    for item in items:
        if item.get_closest_marker("slowest") and not config.getoption(slowest_opion):
            item.add_marker(skip_slowest)
        elif item.get_closest_marker("slow") and not config.getoption(slow_opion):
            item.add_marker(skip_slow)
