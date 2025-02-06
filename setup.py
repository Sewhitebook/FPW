from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
import os
import sys

# Get version from the environment variable passed by GitHub Actions
version = os.getenv('PACKAGE_VERSION', '0.0.0')  # Default to 0.0.0 if not set

extra_compile_args = ["-O3", "-DNDEBUG", "-std=c++17"]

extensions = [
    Extension(
        "fpw",
        ["fpw.pyx", "src/fpw.cpp"],  # Keep the relative path to fpw.cpp
        include_dirs=[numpy.get_include(), 'src'],
        language="c++",
        extra_compile_args=extra_compile_args,
    )
]

setup(
    name="fpw",
    version=version,  # Use the dynamic version here
    ext_modules=cythonize(extensions),
    install_requires=[
        'numpy',
    ],
)

