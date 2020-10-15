#include "linalg_operations.h"

double linalg::sqr(int N, const double* __restrict x) {
    double r = 0;
    for (int i = 0; i < N; ++i)
        r += x[i] * x[i];
    return r;
}

double linalg::dot(int N, const double* __restrict x, const double* __restrict y) {
    double r = 0;
    for (int i = 0; i < N; ++i)
        r += x[i] * y[i];
    return r;
}

void linalg::multiply(int N, double a, double* __restrict x) {
    for (int i = 0; i < N; ++i)
        x[i] = a*x[i];
}

void linalg::axpby(int N, double a, double* __restrict x, double b, const double* __restrict y) {
    for (int i = 0; i < N; ++i)
        x[i] = a*x[i] + b*y[i];
}

void linalg::SpMV(
    int i_begin, int i_end, double* __restrict z,
    const int* __restrict AI, const int* __restrict AJ,
    const double* __restrict A, const double* __restrict x
) {
    for (int i = i_begin; i < i_end; ++i) {
        double r = 0;
        for (int j = AI[i]; j < AI[i+1]; ++j)
            r += A[j] * x[AJ[j]];
        z[i] = r;
    }
}

void linalg::prod(int N, double* __restrict z, const double* __restrict x, const double* __restrict y) {
    for (int i = 0; i < N; ++i)
        z[i] = x[i] * y[i];
}
