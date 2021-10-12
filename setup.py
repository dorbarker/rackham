from setuptools import setup, find_packages
from rackham import __version__

setup(
    name="rackham",
    version=__version__,
    packages=find_packages(),
    install_requires=["biopython", "pandas", "numpy"],
    author="Dillon O.R. Barker",
    author_email="dillon.barker@phac-aspc.gc.ca",
    entry_points={"console_scripts": ["rackham=rackham.rackham:main"]},
)
