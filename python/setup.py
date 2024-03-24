from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        "fpw_cython",
        ["fpw_cython.pyx", "../src/fpw.cpp"],  # Keep the relative path to fpw.cpp
        include_dirs=[numpy.get_include(), '../eigen', '../src'],
        language="c++",
        extra_compile_args=["-O3", "-DNDEBUG"],
    )
]

setup(
    name="fpw",
    ext_modules=cythonize(extensions),
    install_requires=[
        'numpy',
    ],
)