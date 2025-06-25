import io, re
from setuptools import setup

# pull version from the source file
ver = re.search(
    r"^__version__ = ['\"]([^'\"]*)['\"]",
    io.open("PhyREC/_version.py", encoding="utf-8").read(),
    re.M,
).group(1)

if __name__ == "__main__":
    setup(version=ver)
