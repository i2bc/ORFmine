from pathlib import Path

# root path of the project
ROOT_PATH = Path(__file__).resolve().parent.parent

# docker infos
DOCKER_REPOSITORY = "lopesi2bc"
IMAGE_NAME = "orfmine"
IMAGE_TAG = "latest"
DOCKER_IMAGE = f"{DOCKER_REPOSITORY}/{IMAGE_NAME}:{IMAGE_TAG}"

