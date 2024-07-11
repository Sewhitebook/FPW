#include "fpw.h"

/**
 * @brief Given N phase bins, makeIndices maps the bin index that each timestamp falls in.
 * 
 * This function calculates the phase for each timestamp and maps it to a bin index.
 * 
 * @param times Array of timestamps.
 * @param f Single frequency to be tested.
 * @param N_bins Number of phase bins.
 * @param N_dat Number of data points.
 * @return Pointer to an array of bin indices corresponding to each timestamp.
 */
int* makeIndices(double* times, const double f, const int N_bins, const int N_dat){
    int* indices = new int[N_dat];

    for(int i = 0; i < N_dat; i++){
        double phase = times[i] * f;
        phase -= static_cast<int>(phase); // Gets the fractional part of phase
        indices[i] = static_cast<int>(phase * N_bins); // The index for any given phase is the integer part of phase * N_bins
    }
    return indices; // Indices has the length of the timeseries
}

/**
 * @brief Given the indices of the timestamps, deltaChi2 computes the FPW statistic.
 * 
 * This function calculates the FPW statistic by accumulating the inverse variance
 * and the product of the inverse variance and the data for each bin.
 * 
 * @param y Array of observed values.
 * @param ivar Array of inverse variances.
 * @param ind Array of bin indices for each timestamp.
 * @param N_bins Number of phase bins.
 * @param N_dat Number of data points.
 * @return Computed FPW statistic.
 */
double deltaChi2(double* y, double* ivar, int* ind, int N_bins, int N_dat){
    double* VtCinvV = new double[N_bins]();
    double* ytCinvV = new double[N_bins]();
    double* ivar_y = new double[N_dat];

    for (int i = 0; i < N_dat; ++i){
        ivar_y[i] = ivar[i] * y[i];
    }

    double deltaChi = 0.0;
    for (int i = 0; i < N_dat; ++i){
        int index = ind[i];
        VtCinvV[index] += ivar[i]; // Accumulate the inverse variance for each bin
        ytCinvV[index] += ivar_y[i]; // Accumulate the product of the inverse variance and the data for each bin
    }

    for (int i = 0; i < N_bins; ++i){
        double Minv = 1.0 / (2 * VtCinvV[i]); //build the denominator of S_FPW out of the accumulated values
        deltaChi += ytCinvV[i] * Minv * ytCinvV[i]; // Compute the outer sum of S_FPW
    }

    delete[] VtCinvV;
    delete[] ytCinvV;
    delete[] ivar_y;

    return deltaChi;
}

/**
 * @brief Computes the relevant statistic for every test frequency.
 * 
 * This function calculates the FPW statistic for each test frequency
 * by computing the inverse variance and looping over all frequencies.
 * 
 * @param t Vector of time points.
 * @param y Vector of observed values.
 * @param dy Vector of uncertainties in the observed values.
 * @param freqs Vector of test frequencies.
 * @param N_bins Number of bins to use in the computation.
 * @return Vector of computed statistics for each frequency.
 */
vector<double> runFPW(vector<double>& t, vector<double>& y, vector<double>& dy, const vector<double>& freqs, int N_bins){
    int N_freqs = freqs.size();
    int N_dat = t.size();

    // Compute the inverse variance
    vector<double> ivar(N_dat);
    for (int i = 0; i < N_dat; ++i){
        ivar[i] = 1.0 / (dy[i] * dy[i]);
    }

    // Loop over all frequencies and compute the FPW statistic
    vector<double> deltaChiArr(N_freqs);
    for (int i = 0; i < N_freqs; ++i){
        int* indices = makeIndices(t.data(), freqs[i], N_bins, N_dat);
        deltaChiArr[i] = deltaChi2(y.data(), ivar.data(), indices, N_bins, N_dat);
        delete[] indices;
    }

    return deltaChiArr;
}

/**
 * @brief Computes the relevant statistic for multiple light curves sharing the same timestamps.
 * 
 * This function calculates the FPW statistic for each test frequency
 * for multiple light curves by computing the inverse variance
 * and looping over all frequencies.
 * 
 * @param t Vector of time points.
 * @param y Vector of vectors of observed values for multiple light curves.
 * @param dy Vector of vectors of uncertainties in the observed values for multiple light curves.
 * @param freqs Vector of test frequencies.
 * @param N_bins Number of bins to use in the computation.
 * @return Vector of vectors of computed statistics for each frequency for each light curve.
 */
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