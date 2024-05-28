#ifndef FPW_H
#define FPW_H

#include <vector>
using std::vector;

int* makeIndices(double* times, const double f, const int N_bins, const int size);
double deltaChi2(double* y, double* ivar, int* ind, int N_bins, int N_dat);
vector<double> runFPW(vector<double>& t, vector<double>& y, vector<double>& dy, const vector<double>& freqs, int N_bins);
// vector<double> testFPW(const vector<double>& t, const vector<double>& y, const vector<double>& ivar, const vector<double>& freqs, int N_bins);

#endif // FPW_H