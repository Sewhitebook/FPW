#include "fpw.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::VectorXi;

VectorXd mapArrayToVectorXd(double* data, int size) {
    return Eigen::Map<Eigen::VectorXd>(data, size);
}

VectorXi makeIndices(const VectorXd& times, const double f, const int N_bins){
    VectorXd tf = times * f;
    tf = tf.array() - tf.array().cast<int>().cast<double>();
    VectorXi indices = (tf * N_bins).cast<int>();
    return indices;
}

double deltaChi2(const VectorXd& y, const VectorXd& ivar, const VectorXi& ind, int N_bins, int N_dat){
    VectorXd acc = VectorXd::Zero(N_bins);
    VectorXd ytCinvV = VectorXd::Zero(N_bins);
    VectorXd ivar_y = ivar.array() * y.array();

    double deltaChi = 0.0;
    for (int i = 0; i < N_dat; ++i){
        int index = ind[i];
        acc[index] += ivar[i];
        ytCinvV[index] += ivar_y[i];
    }

    for (int i = 0; i < N_bins; ++i){
        double Minv = 1.0 / (acc[i] + 1.0);
        deltaChi += ytCinvV[i] * Minv * ytCinvV[i];
    }
    
    return deltaChi;
}

VectorXd runFPW(const VectorXd& t, const VectorXd& y, const VectorXd& ivar, const VectorXd& freqs, int N_bins){
    int N_freqs = freqs.size();
    int N_dat = t.size();

    VectorXd deltaChiArr = VectorXd::Zero(N_freqs);
    for (int i = 0; i < N_freqs; ++i){
        VectorXi indices = makeIndices(t, freqs[i], N_bins);
        deltaChiArr[i] = deltaChi2(y, ivar, indices, N_bins, N_dat);
    }

    return deltaChiArr;
}