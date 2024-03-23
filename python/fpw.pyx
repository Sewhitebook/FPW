# distutils: language = c++
# distutils: include_dirs = ../eigen

cimport numpy as np
import numpy as np

cdef extern from "fpw.h":
    cdef cppclass VectorXd "Eigen::VectorXd":
        VectorXd()
        VectorXd(int size)
        double& operator[](int i)
        int size()

    cdef cppclass VectorXi "Eigen::VectorXi":
        VectorXi()
        VectorXi(int size)
        int& operator[](int i)
        int size()

    VectorXd mapArrayToVectorXd(double* data, int size)
    VectorXd runFPW(VectorXd t, VectorXd y, VectorXd ivar, VectorXd freqs, int N_bins)

def py_runFPW(np.ndarray[double, ndim=1] t, np.ndarray[double, ndim=1] y, np.ndarray[double, ndim=1] ivar, np.ndarray[double, ndim=1] freqs, int N_bins):
    cdef:
        VectorXd c_t = mapArrayToVectorXd(&t[0], t.size)
        VectorXd c_y = mapArrayToVectorXd(&y[0], y.size)
        VectorXd c_ivar = mapArrayToVectorXd(&ivar[0], ivar.size)
        VectorXd c_freqs = mapArrayToVectorXd(&freqs[0], freqs.size)
        VectorXd c_result

    c_result = runFPW(c_t, c_y, c_ivar, c_freqs, N_bins)

    result = np.zeros(c_result.size())
    for i in range(c_result.size()):
        result[i] = c_result[i]

    return result