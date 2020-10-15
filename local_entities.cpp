#include "local_entities.h"

Matrix& Matrix::fill(
    const Space &sp, Portrait &&pt,
    function<double(int,int)> a_ij, function<double(int,double)> a_ii
) {
    AI = std::move(pt.AI);
    AJ = std::move(pt.AJ);
    A.resize(AJ.size());
    int N = AI.size()-1;

    for (int i = 0; i < N; ++i) {
        double s = 0;
        int k = -1;
        int i_gl = sp.l2g[i];

        for (int j = AI[i]; j < AI[i+1]; ++j)
            if (AJ[j] != i) {
                int j_gl = sp.l2g[AJ[j]];
                double r = a_ij(i_gl, j_gl);
                A[j] = r;
                s += fabs(r);
            } else
                k = j;

        if (k != -1)
            A[k] = a_ii(i_gl, s);
    }

    return *this;
}


vector<double> fill_vector(
    const Space& sp, function<double(int)> b_i
) {
    int N = sp.N_own;
    vector<double> b(N);

    for (int i = 0; i < N; ++i) {
        int i_gl = sp.l2g[i];
        b[i] = b_i(i_gl);
    }

    return b;
}


vector<double> jacobi_preconditioner(const Matrix& A) {
    int N = A.AI.size() - 1;
    vector<double> z(N);

    for (int i = 0; i < N; ++i) {
        double s = -1;
        int jl = A.AI[i];
        int jr = A.AI[i+1]-1;
        while (jl < jr) {
            int jm = (jl+jr)/2;
            if (i <= A.AJ[jm])
                jr = jm;
            else
                jl = jm+1;
        }

        if (i == A.AJ[jl])
            s = A.A[jl];

        if (s > 0)
            z[i] = 1/s;
        else
            { /* error */ }
    }
    return z;
}
