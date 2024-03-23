#ifndef FPW_H
#define FPW_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::VectorXi;

VectorXi makeIndices(const VectorXd& times, const double f, const int N_bins);
double deltaChi2(const VectorXd& y, const VectorXd& ivar, const VectorXi& ind, int N_bins, int N_dat);
VectorXd runFPW(const VectorXd& t, const VectorXd& y, const VectorXd& ivar, const VectorXd& freqs, int N_bins);
VectorXd mapArrayToVectorXd(double* data, int size);

#endif // FPW_H