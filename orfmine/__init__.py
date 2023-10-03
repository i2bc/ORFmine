from pathlib import Path
import os

# root path of the project
ROOT_PATH = Path(__file__).resolve().parent.parent


IMAGE_TAGS = ["latest", "2.0"]

# docker infos
DOCKER_REPOSITORY = os.getenv("DOCKER_REPOSITORY", "lopesi2bc")
IMAGE_NAME = os.getenv("IMAGE_NAME", "orfmine")
IMAGE_TAG = os.getenv("IMAGE_TAG", IMAGE_TAGS[-1])

DOCKER_IMAGE = os.getenv("DOCKER_IMAGE", f"{DOCKER_REPOSITORY}/{IMAGE_NAME}:{IMAGE_TAG}")

