from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension(
        "fpw",
        ["fpw.pyx", "../src/fpw.cpp"],
        include_dirs=[numpy.get_include(), '../eigen', '../src'],
        language="c++"
    )
]

setup(
    name="fpw",
    ext_modules=cythonize(extensions),
    install_requires=[
        'numpy',
    ],
)