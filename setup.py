#!/usr/bin/env python3

from setuptools import setup
import os


# We store our version number in a simple text file:
version_path = os.path.join(
    os.path.dirname(__file__), "chemical_composition", "version.txt"
)
with open(version_path, "r") as version_file:
    chemical_composition_version = version_file.read().strip()

setup(
    name="chemical_composition",
    version=chemical_composition_version,
    packages=["chemical_composition"],
    package_dir={"chemical_composition": "chemical_composition"},
    description="chemical_composition",
    package_data={
        "chemical_composition": [
            "version.txt",
        ]
    },
    long_description="Unimod Mapper for Proteomics tools",
    author="... and Christian Fufezan",
    author_email="christian@fufezan.net",
    url="https://github.com/computational-ms/unimod-mapper",
    license="Lesser GNU General Public License (LGPL)",
    platforms="any that supports python 3.6",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: SunOS/Solaris",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Education",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
)
