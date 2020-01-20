import pytest


def pytest_addoption(parser):
    parser.addoption(
            "--runslow", action="store_true", default=False, help="run slow tests"
            )
    parser.addoption(
            "--run-very-slow", action="store_true", default=False, help="run very slow tests"
            )


    def pytest_configure(config):
        config.addinivalue_line("markers", "slow: mark test as slow to run")
        config.addinivalue_line("markers", "veryslow: mark test as very slow to run")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--runslow"):
        skip_slow = pytest.mark.skip(reason="need --runslow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)

    if not config.getoption("--run-very-slow"):
        skip_veryslow = pytest.mark.skip(reason="need --run-very-slow option to run")
        for item in items:
            if "veryslow" in item.keywords:
                item.add_marker(skip_veryslow)

