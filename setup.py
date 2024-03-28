from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
import sys

extra_compile_args = ["-O3", "-DNDEBUG"]
if sys.platform == "darwin":  # macOS
    extra_compile_args += ["-std=c++11", "-Wno-nullability-completeness"]

extensions = [
    Extension(
        "fpw",
        ["fpw.pyx", "src/fpw.cpp"],  # Keep the relative path to fpw.cpp
        include_dirs=[numpy.get_include(), 'include/eigen', 'src'],
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