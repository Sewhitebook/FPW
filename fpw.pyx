# distutils: language = c++
# distutils: sources = src/fpw.cpp

import numpy as np
cimport numpy as cnp
from libcpp.vector cimport vector

cdef extern from "fpw.h":
    vector[double] runFPW(vector[double]& t, vector[double]& y, vector[double]& dy, vector[double]& freqs, int N_bins)
    vector[vector[double]] runFPWMulti(vector[double]& t, vector[vector[double]]& y, vector[vector[double]]& dy, vector[double]& freqs, int N_bins)

def run_fpw(cnp.ndarray[double, ndim=1, mode="c"] t, cnp.ndarray[double, ndim=1, mode="c"] y, cnp.ndarray[double, ndim=1, mode="c"] dy, cnp.ndarray[double, ndim=1, mode="c"] freqs, int N_bins):
    cdef vector[double] t_vec = <vector[double]&>t
    cdef vector[double] y_vec = <vector[double]&>y
    cdef vector[double] dy_vec = <vector[double]&>dy
    cdef vector[double] freqs_vec = <vector[double]&>freqs

    cdef vector[double] result = runFPW(t_vec, y_vec, dy_vec, freqs_vec, N_bins)

    # Convert the C++ vector to a numpy array
    return np.array(result)

def run_fpw_multi(cnp.ndarray[double, ndim=1, mode="c"] t, cnp.ndarray[double, ndim=2, mode="c"] y, cnp.ndarray[double, ndim=2, mode="c"] dy, cnp.ndarray[double, ndim=1, mode="c"] freqs, int N_bins):
    cdef vector[double] t_vec = <vector[double]&>t

    assert y.shape[1] == t.shape[0], "The number of columns in y must match the length of t"
    assert dy.shape[1] == t.shape[0], "The number of columns in dy must match the length of t"

    cdef vector[vector[double]] y_vec = <vector[vector[double]]&>y
    cdef vector[vector[double]] dy_vec = <vector[vector[double]]&>dy
    cdef vector[double] freqs_vec = <vector[double]&>freqs

    cdef vector[vector[double]] result = runFPWMulti(t_vec, y_vec, dy_vec, freqs_vec, N_bins)

    # Convert the C++ vector of vectors to a numpy array
    return np.array([[result[i][j] for j in range(len(result[i]))] for i in range(len(result))])