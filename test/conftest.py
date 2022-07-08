from email.policy import default

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--skip-slow",
        action="store_true",
        default=False,
        help="skip slow tests",
    )
    parser.addoption(
        "--parallel",
        type=int,
        default=1,
        help="run tests in parallel processes",
    )

@pytest.fixture
def the_number_of_process(request):
    return request.config.getoption("--parallel")
