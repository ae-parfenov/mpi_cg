#include "shared_operations.h"

double shared::min(const Space& sp, Gate& gt, const vector<double>& x) {
    int N = sp.N_own;
    double r = x[0];
    for (int i = 1; i < N; ++i)
        if (x[i] < r) r = x[i];
    r = gt.allreduce(r, MPI_MIN);
    return r;
}

double shared::max(const Space& sp, Gate& gt, const vector<double>& x) {
    int N = sp.N_own;
    double r = x[0];
    for (int i = 1; i < N; ++i)
        if (x[i] > r) r = x[i];
    r = gt.allreduce(r, MPI_MAX);
    return r;
}

double shared::dot(
    const Space& sp, Gate& gt,
    const vector<double>& x, const vector<double>& y, double* t
) {
    int N = sp.N_own;
    double r = 0;
    double t0 = MPI_Wtime();

    if (&x == &y)
        r = linalg::sqr(N, x.data());
    else
        r = linalg::dot(N, x.data(), y.data());
    r = gt.allreduce(r, MPI_SUM);

    if(t) *t += MPI_Wtime() - t0;
    return r;
}

void shared::axpby(
    const Space& sp, Gate& /*gt*/,
    double a, vector<double>& x, double b, const vector<double>& y, double* t
) {
    int N = sp.N_own;
    double t0 = MPI_Wtime();

    if (&x == &y)
        linalg::multiply(N, a+b, x.data());
    else
        linalg::axpby(N, a, x.data(), b, y.data());

    if(t) *t += MPI_Wtime() - t0;
}

void shared::SpMV(
        const Space& sp, Gate& gt,
        vector<double>& z, const Matrix& A, vector<double>& x, double* t
) {
    double t0 = MPI_Wtime();
    gt.exchange(x);
    linalg::SpMV(0, sp.N_inner, z.data(), A.AI.data(), A.AJ.data(), A.A.data(), x.data());
    gt.waitall();
    linalg::SpMV(sp.N_inner, sp.N_own, z.data(), A.AI.data(), A.AJ.data(), A.A.data(), x.data());

    if(t) *t += MPI_Wtime() - t0;
}

void shared::prod(
    const Space& sp, Gate& /*gt*/,
    vector<double>& z, const vector<double>& x, const vector<double>& y, double* t
) {
    double t0 = MPI_Wtime();

    int N = sp.N_own;
    linalg::prod(N, z.data(), x.data(), y.data());

    if(t) *t += MPI_Wtime() - t0;
}


vector<double> shared::solve_cg(
    const Space& sp, Gate& gt, const Matrix& A, const vector<double>& b, double tol,
    bool print_log, Stats* st
) {
    double t = MPI_Wtime();
    double *t_dot, *t_axpby, *t_spmv, *t_prod;

    if (st) {
        st->times.resize(4);
        st->calls.resize(4);
        t_dot = &(st->times[0]);
        t_axpby = &(st->times[1]);
        t_spmv = &(st->times[2]);
        t_prod = &(st->times[3]);
        *t_dot = *t_axpby = *t_spmv = *t_prod = 0;
    } else
        t_dot = t_axpby = t_spmv = t_prod = nullptr;

    vector<double> B = jacobi_preconditioner(A);
    vector<double> x(sp.N_own, 0);
    vector<double> p(sp.N_local, 0);
    vector<double> r(b);
    vector<double> z(sp.N_own);
    vector<double> q(sp.N_own);

    double rho, alpha, beta;
    double l2, eps;
    int k = 0;

    l2 = sqrt(dot(sp, gt, b, b, t_dot));
    eps = tol * l2;
    rho = 1;

    if (print_log && gt.my_rank == 0) {
        printf("\n-- CG Start -- |Ax-b|\n");
        printf(" iteration %d : %f \n", k, l2);
        fflush(stdout);
    }

    while (l2 > eps) {
        ++k;
        prod(sp, gt, z, B, r, t_prod);
        beta = rho;
        rho = dot(sp, gt, r, z, t_dot);
        axpby(sp, gt, rho/beta, p, 1, z, t_axpby);
        SpMV(sp, gt, q, A, p, t_spmv);
        alpha = rho / dot(sp, gt, p, q, t_dot);
        axpby(sp, gt, 1, x, alpha, p, t_axpby);
        axpby(sp, gt, 1, r, -alpha, q, t_axpby);
        l2 = sqrt(dot(sp, gt, r, r, t_dot));

        if (print_log && gt.my_rank == 0) {
            printf(" iteration %d : %f \n", k, l2);
            fflush(stdout);
        }
    }

    t = MPI_Wtime() - t;

    if (print_log && gt.my_rank == 0) {
        printf("final |Ax-b|/|b|: %g\n", l2 / (eps/tol));
        printf("total time : %f s\n", t);
        fflush(stdout);
    }

    if (st) {
        st->total_time = t;
        st->calls[0] = 2 + k*3;
        st->calls[1] = 1 + k*3;
        st->calls[2] = 1 + k;
        st->calls[3] = k > 0 ? k : 1;
    }

    return x;
}
