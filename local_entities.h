#ifndef LOCAL_H
#define LOCAL_H

#include "spaces_connection.h"

#include <math.h>

#include <vector>
#include <functional>

using std::vector;
using std::function;

struct Matrix {
    vector<int> AI, AJ;
    vector<double> A;

    Matrix& fill(
        const Space &sp,
        Portrait &&pt,
        function<double(int,int)> a_ij = [](int i, int j)->double { return cos(1.0*i*j); },
        function<double(int,double)> a_ii = [](int /*i*/, double s)->double { return 1.5*s; }
    );
};


vector<double> fill_vector(
    const Space& sp,
    function<double(int)> b_i = [](int i)->double { return cos(1.0*i*i); }
);


vector<double> jacobi_preconditioner(const Matrix& A);

#endif // LOCAL_H
