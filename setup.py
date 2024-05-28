from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
import sys

extra_compile_args = ["-O3", "-DNDEBUG", "-std=c++17"]

extensions = [
    Extension(
        "fpw",
        ["fpw.pyx", "src/fpw.cpp"],  # Keep the relative path to fpw.cpp
        # include_dirs=[numpy.get_include(), 'include/eigen', 'src'],
        include_dirs=[numpy.get_include(), 'src'],
        language="c++",
        extra_compile_args=extra_compile_args,
    )
]

setup(
    name="fpw",
    ext_modules=cythonize(extensions),
    install_requires=[
        'numpy',
    ],
)
