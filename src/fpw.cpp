#include "fpw.h"

using Eigen::VectorXd;
using Eigen::VectorXi;
using namespace::std;

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

VectorXd runFPW(const VectorXd& y, const VectorXd& ivar, const VectorXd& freqs, const VectorXd& t, int N_bins){
    int N_freqs = freqs.size();
    int N_dat = t.size();

    VectorXd deltaChiArr = VectorXd::Zero(N_freqs);
    for (int i = 0; i < N_freqs; ++i){
        VectorXi indices = makeIndices(t, freqs[i], N_bins);
        deltaChiArr[i] = deltaChi2(y, ivar, indices, N_bins, N_dat);
    }

    return deltaChiArr;
}

// int main(){
//     vector<double> t, y, dy, freqs;
//     string line, value;

//     // Read t, y, dy from the CSV file
//     ifstream tydyFile("/data/sams_datasets/tydy.csv");
//     while (getline(tydyFile, line)) {
//         stringstream s(line);
//         getline(s, value, ',');
//         t.push_back(stod(value));
//         getline(s, value, ',');
//         y.push_back(stod(value));
//         getline(s, value, ',');
//         dy.push_back(stod(value));
//     }

//     Eigen::VectorXd t_eigen = Eigen::VectorXd::Map(t.data(), t.size());
//     Eigen::VectorXd y_eigen = Eigen::VectorXd::Map(y.data(), y.size());
//     Eigen::VectorXd dy_eigen = Eigen::VectorXd::Map(dy.data(), dy.size());

//     // read frequencies from fgrid.csv
//     ifstream freqFile("/data/sams_datasets/fgrid.csv");
//     while (getline(freqFile, line)) {
//         freqs.push_back(stod(line));
//     }

//     // cast freqs to Eigen::VectorXd
//     Eigen::VectorXd freqs_eigen = Eigen::VectorXd::Map(freqs.data(), freqs.size());

//     int N_freqs = freqs.size();
//     int N_dat = t.size();
//     int N_bins = 10;

//     // define ivar by dy
//     Eigen::VectorXd ivar = 1.0 / (dy_eigen.array() * dy_eigen.array());

//     // loop each frequency and find the time to compute a whole grid
//     Eigen::VectorXd deltaChiArr = runFPW(y_eigen, ivar, freqs_eigen, t_eigen, N_bins);

//     std::ofstream outputFile("/data/sams_datasets/FWP_out.csv");
//     for (int i = 0; i < N_freqs; ++i){
//         outputFile << freqs_eigen[i] << "," << deltaChiArr[i] << "\n";
//     }
//     outputFile.close();
//     return 0;
// }