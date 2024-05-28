# distutils: language = c++
# distutils: sources = src/fpw.cpp

import numpy as np
cimport numpy as cnp
from libcpp.vector cimport vector

cdef extern from "fpw.h":
    vector[double] runFPW(vector[double]& t, vector[double]& y, vector[double]& dy, vector[double]& freqs, int N_bins)

def run_fpw(cnp.ndarray[double, ndim=1, mode="c"] t, cnp.ndarray[double, ndim=1, mode="c"] y, cnp.ndarray[double, ndim=1, mode="c"] dy, cnp.ndarray[double, ndim=1, mode="c"] freqs, int N_bins):
    cdef vector[double] t_vec = <vector[double]&>t
    cdef vector[double] y_vec = <vector[double]&>y
    cdef vector[double] dy_vec = <vector[double]&>dy
    cdef vector[double] freqs_vec = <vector[double]&>freqs

    cdef vector[double] result = runFPW(t_vec, y_vec, dy_vec, freqs_vec, N_bins)

    # Convert the C++ vector to a numpy array
    return np.array(result)