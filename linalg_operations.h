#ifndef LINALG_H
#define LINALG_H

namespace linalg {

    double sqr(int N, const double* __restrict x);

    double dot(int N, const double* __restrict x, const double* __restrict y);

    void multiply(int N, double a, double* __restrict x);

    void axpby(int N, double a, double* __restrict x, double b, const double* __restrict y);

    void SpMV(
        int i_begin, int i_end, double* __restrict z,
        const int* __restrict AI, const int* __restrict AJ,
        const double* __restrict A, const double* __restrict x
    );

    void prod(int N, double* __restrict z, const double* __restrict x, const double* __restrict y);
}

#endif // LINALG_H
