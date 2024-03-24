#include "fpw.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using std::vector;

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

vector<double> runFPW(const vector<double>& t, const vector<double>& y, const vector<double>& ivar, const vector<double>& freqs, int N_bins){
    VectorXd t_eigen = VectorXd::Map(t.data(), t.size());
    VectorXd y_eigen = VectorXd::Map(y.data(), y.size());
    VectorXd ivar_eigen = VectorXd::Map(ivar.data(), ivar.size());
    VectorXd freqs_eigen = VectorXd::Map(freqs.data(), freqs.size());

    int N_freqs = freqs_eigen.size();
    int N_dat = t_eigen.size();

    vector<double> deltaChiArr(N_freqs);
    for (int i = 0; i < N_freqs; ++i){
        VectorXi indices = makeIndices(t_eigen, freqs_eigen[i], N_bins);
        deltaChiArr[i] = deltaChi2(y_eigen, ivar_eigen, indices, N_bins, N_dat);
    }

    return deltaChiArr;
}