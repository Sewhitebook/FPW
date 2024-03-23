# distutils: language = c++
# distutils: include_dirs = ../eigen

cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free

cdef extern from "fpw.h":
    cdef cppclass VectorXd "Eigen::VectorXd":
        VectorXd()
        VectorXd(int size)
        double operator[](int i)
        int size()

    cdef cppclass VectorXi "Eigen::VectorXi":
        VectorXi()
        VectorXi(int size)
        int operator[](int i)
        int size()

    VectorXi makeIndices(VectorXd times, double f, int N_bins)
    double deltaChi2(VectorXd y, VectorXd ivar, VectorXi ind, int N_bins, int N_dat)
    VectorXd runFPW(VectorXd y, VectorXd ivar, VectorXd freqs, VectorXd t, int N_bins)
    VectorXd mapArrayToVectorXd(double* data, int size)

def py_runFPW(np.ndarray[double, ndim=1] y, np.ndarray[double, ndim=1] ivar, np.ndarray[double, ndim=1] freqs, np.ndarray[double, ndim=1] t, int N_bins):
    cdef:
        VectorXd c_y = VectorXd(y.size)
        VectorXd c_ivar = VectorXd(ivar.size)
        VectorXd c_freqs = VectorXd(freqs.size)
        VectorXd c_t = VectorXd(t.size)
        VectorXd c_result

    # Copy data from numpy arrays to VectorXd objects
    cdef double* y_data = <double*> malloc(y.size * sizeof(double))
    cdef double* ivar_data = <double*> malloc(ivar.size * sizeof(double))
    cdef double* freqs_data = <double*> malloc(freqs.size * sizeof(double))
    cdef double* t_data = <double*> malloc(t.size * sizeof(double))

    for i in range(y.size):
        y_data[i] = y[i]
        ivar_data[i] = ivar[i]
        freqs_data[i] = freqs[i]
        t_data[i] = t[i]

    c_y = mapArrayToVectorXd(y_data, y.size)
    c_ivar = mapArrayToVectorXd(ivar_data, ivar.size)
    c_freqs = mapArrayToVectorXd(freqs_data, freqs.size)
    c_t = mapArrayToVectorXd(t_data, t.size)

    c_result = runFPW(c_y, c_ivar, c_freqs, c_t, N_bins)

    result = np.zeros(c_result.size())
    for i in range(c_result.size()):
        result[i] = c_result[i]

    # Free the allocated memory
    free(y_data)
    free(ivar_data)
    free(freqs_data)
    free(t_data)

    return result