#ifndef SHARED_H
#define SHARED_H

#include "spaces_connection.h"
#include "local_entities.h"
#include "linalg_operations.h"

#include <stdio.h>
#include <mpi.h>

#include<vector>

using std::vector;

struct Stats {
    vector<int> calls;
    vector<double> times;
    double total_time;

    Stats() : calls(4,0), times(4,0.0), total_time(0) {    }
};


namespace shared {

    double min(const Space& sp, Gate& gt, const vector<double>& x);

    double max(const Space& sp, Gate& gt, const vector<double>& x);

    double dot(
        const Space& sp, Gate& gt,
        const vector<double>& x, const vector<double>& y, double* t = nullptr
    );

    void axpby(
        const Space& sp, Gate& gt,
        double a, vector<double>& x, double b, const vector<double>& y, double* t = nullptr
    );

    void SpMV(
        const Space& sp, Gate& gt,
        vector<double>& z, const Matrix& A, vector<double>& x, double* t = nullptr
    );

    void prod(
        const Space& sp, Gate& gt,
        vector<double>& z, const vector<double>& x, const vector<double>& y, double* t = nullptr
    );


    vector<double> solve_cg(
        const Space& sp, Gate& gt, const Matrix& A, const vector<double>& b, double tol,
        bool print_log = false, Stats* st = nullptr
    );
}

#endif // SHARED_H
