#include "fpw.h"

/*Given N phase bins, makeIndices maps the bin index that each timestamp falls in*/
int* makeIndices(double* times, const double f, const int N_bins, const int size){
    int* indices = new int[size];

    for(int i = 0; i < size; i++){
        double tf = times[i] * f;
        tf -= static_cast<int>(tf);
        indices[i] = static_cast<int>(tf * N_bins);
    }
    return indices;
}

double deltaChi2(double* y, double* ivar, int* ind, int N_bins, int N_dat){
    double* acc = new double[N_bins]();
    double* ytCinvV = new double[N_bins]();
    double* ivar_y = new double[N_dat];

    for (int i = 0; i < N_dat; ++i){
        ivar_y[i] = ivar[i] * y[i];
    }

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

    delete[] acc;
    delete[] ytCinvV;
    delete[] ivar_y;

    return deltaChi;
}

vector<double> runFPW(vector<double>& t, vector<double>& y, vector<double>& dy, const vector<double>& freqs, int N_bins){
    int N_freqs = freqs.size();
    int N_dat = t.size();

    vector<double> ivar(N_dat);
    for (int i = 0; i < N_dat; ++i){
        ivar[i] = 1.0 / (dy[i] * dy[i]);
    }

    vector<double> deltaChiArr(N_freqs);
    for (int i = 0; i < N_freqs; ++i){
        int* indices = makeIndices(t.data(), freqs[i], N_bins, N_dat);
        deltaChiArr[i] = deltaChi2(y.data(), ivar.data(), indices, N_bins, N_dat);
        delete[] indices;
    }

    return deltaChiArr;
}

vector<vector<double>> runFPWMulti(vector<double>& t, vector<vector<double>>& y, vector<vector<double>>& dy, const vector<double>& freqs, int N_bins){
    int N_freqs = freqs.size();
    int N_dat = t.size();
    int N_curves = y.size();

    vector<vector<double>> ivar(N_curves, vector<double>(N_dat));
    for (int j = 0; j < N_curves; ++j){
        for (int i = 0; i < N_dat; ++i){
            ivar[j][i] = 1.0 / (dy[j][i] * dy[j][i]);
        }
    }

    vector<vector<double>> deltaChiArr(N_curves, vector<double>(N_freqs));
    for (int i = 0; i < N_freqs; ++i){
        int* indices = makeIndices(t.data(), freqs[i], N_bins, N_dat);
        for (int j = 0; j < N_curves; ++j){
            deltaChiArr[j][i] = deltaChi2(y[j].data(), ivar[j].data(), indices, N_bins, N_dat);
        }
        delete[] indices;
    }
    return deltaChiArr;
}