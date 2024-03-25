#ifndef FPW_H
#define FPW_H

#include <vector>
#include "Eigen/Dense"

using Eigen::VectorXd;
using Eigen::VectorXi;
using std::vector;

VectorXi makeIndices(const VectorXd& times, const double f, const int N_bins);
double deltaChi2(const VectorXd& y, const VectorXd& ivar, const VectorXi& ind, int N_bins, int N_dat);
vector<double> runFPW(const vector<double>& t, const vector<double>& y, const vector<double>& ivar, const vector<double>& freqs, int N_bins);
vector<double> testFPW(const vector<double>& t, const vector<double>& y, const vector<double>& ivar, const vector<double>& freqs, int N_bins);

#endif // FPW_H