import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--plot", action="store_true", default=False
    )

@pytest.fixture
def cmd_plot(request):
    return request.config.getoption("--plot")