from importlib import resources
from pathlib import Path


RESOURCE_DIR = resources.files(__package__)


def get_resource_path(file: str) -> Path:
    return Path(RESOURCE_DIR.joinpath(file)).absolute()
