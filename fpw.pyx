# distutils: language = c++
# distutils: sources = src/fpw.cpp

import numpy as np
cimport numpy as cnp
from libcpp.vector cimport vector

cdef extern from "fpw.h":
    vector[double] runFPW(vector[double]& t, vector[double]& y, vector[double]& dy, vector[double]& freqs, int N_bins)
    vector[vector[double]] runFPWMulti(vector[double]& t, vector[vector[double]]& y, vector[vector[double]]& dy, vector[double]& freqs, int N_bins)
    vector[double] phaseEntropy(vector[double]& t, vector[double]& freqs, int N_bins)

def run_fpw(t, y, dy, freqs, int N_bins):
    if N_bins <= 0:
        raise ValueError("N_bins must be positive")

    t_arr = np.asarray(t, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    dy_arr = np.asarray(dy, dtype=np.float64)
    f_arr = np.asarray(freqs, dtype=np.float64)

    if t_arr.ndim != 1 or y_arr.ndim != 1 or dy_arr.ndim != 1 or f_arr.ndim != 1:
        raise ValueError("t, y, dy, and freqs must be 1D arrays")
    if not (t_arr.shape[0] == y_arr.shape[0] == dy_arr.shape[0]):
        raise ValueError("t, y, and dy must have the same length")

    t_arr = t_arr - np.min(t_arr) #Require t to start at 0
    t = np.ascontiguousarray(t_arr)
    y = np.ascontiguousarray(y_arr)
    dy = np.ascontiguousarray(dy_arr)
    freqs = np.ascontiguousarray(f_arr)

    cdef vector[double] t_vec = <vector[double]&>t
    cdef vector[double] y_vec = <vector[double]&>y
    cdef vector[double] dy_vec = <vector[double]&>dy
    cdef vector[double] freqs_vec = <vector[double]&>freqs
    cdef vector[double] result = runFPW(t_vec, y_vec, dy_vec, freqs_vec, N_bins)
    return np.array(result)

def run_fpw_multi(t, y, dy, freqs, int N_bins):
    if N_bins <= 0:
        raise ValueError("N_bins must be positive")

    t_arr = np.asarray(t, dtype=np.float64)
    y_arr = np.asarray(y, dtype=np.float64)
    dy_arr = np.asarray(dy, dtype=np.float64)
    f_arr = np.asarray(freqs, dtype=np.float64)

    if t_arr.ndim != 1:
        raise ValueError("t must be 1D")
    if y_arr.ndim != 2 or dy_arr.ndim != 2:
        raise ValueError("y and dy must be 2D arrays shaped (n_series, n_points)")
    if y_arr.shape != dy_arr.shape:
        raise ValueError("y and dy must have the same shape")
    if y_arr.shape[1] != t_arr.shape[0]:
        raise ValueError("y and dy second dimension must match length of t")
    if f_arr.ndim != 1:
        raise ValueError("freqs must be 1D")

    t_arr = t_arr - np.min(t_arr) #Require t to start at 0
    t = np.ascontiguousarray(t_arr)
    y = np.ascontiguousarray(y_arr)
    dy = np.ascontiguousarray(dy_arr)
    freqs = np.ascontiguousarray(f_arr)

    cdef vector[double] t_vec = <vector[double]&>t
    cdef vector[vector[double]] y_vec = <vector[vector[double]]&>y
    cdef vector[vector[double]] dy_vec = <vector[vector[double]]&>dy
    cdef vector[double] freqs_vec = <vector[double]&>freqs
    cdef vector[vector[double]] result = runFPWMulti(t_vec, y_vec, dy_vec, freqs_vec, N_bins)
    return np.array([[result[i][j] for j in range(len(result[i]))] for i in range(len(result))])

def phase_entropy(t, freqs, int N_bins):
    if N_bins <= 0:
        raise ValueError("N_bins must be positive")

    t_arr = np.asarray(t, dtype=np.float64)
    f_arr = np.asarray(freqs, dtype=np.float64)

    if t_arr.ndim != 1 or f_arr.ndim != 1:
        raise ValueError("t and freqs must be 1D arrays")

    t_arr = t_arr - np.min(t_arr) #Require t to start at 0
    t = np.ascontiguousarray(t_arr)
    freqs = np.ascontiguousarray(f_arr)

    cdef vector[double] t_vec = <vector[double]&>t
    cdef vector[double] freqs_vec = <vector[double]&>freqs
    cdef vector[double] result = phaseEntropy(t_vec, freqs_vec, N_bins)
    return np.array(result)