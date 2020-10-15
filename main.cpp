#include "global_entities.h"
#include "spaces_connection.h"
#include "local_entities.h"
#include "linalg_operations.h"
#include "shared_operations.h"

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

MPI_Comm MCW = MPI_COMM_WORLD;

#include<vector>

using std::vector;

void print_help() {
    printf("Params: Nx, Ny, Nz, k1, k2, Px, Py [, tol]\n");
    printf(" Nx, Ny, Nz - int > 0 - numbers of layers of nodes along the X,Y,Z axes\n");
    printf(" k1, k2 - int >= 0 - numbers of undivided and divided layers per cycle\n");
    printf(" Px, Py - int > 0, Px * Py = number of processes\n");
    printf(" tol - float > 0, optional (default = 1e-6) - max |Ax-b|/|b|\n");
}


void exchg_report(const Space& sp, const Gate& gt) {
    for (int p = 0; p < gt.num_procs; ++p) {
        if (gt.my_rank == p) {
            printf("\nProcess %d reporting\n", gt.my_rank);
            printf("N| inner : %d, interface : %d, own : %d, halo : %d, local : %d",
                    sp.N_inner, sp.N_interface, sp.N_own, sp.N_halo, sp.N_local);

            for (int i = 0; i < sp.N_local; ++i) {
                if (i == 0) {
                    printf("\nMy inner nodes");
                    if (sp.N_inner > 0) printf(" (0 .. %d)", sp.N_inner-1);
                    printf(":");
                }
                if (i == sp.N_inner) {
                    printf("\nMy interface nodes");
                    if (sp.N_interface > 0) printf(" (%d .. %d)", sp.N_inner, sp.N_own-1);
                    printf(":");
                }
                if (i == sp.N_own) {
                    printf("\nMy halo nodes");
                    if (sp.N_halo > 0) printf(" (%d .. %d)", sp.N_own, sp.N_local-1);
                    printf(":");
                }
                int gl = sp.l2g[i];
                printf(" %d", gl);
                if (i != sp.g2l.at(gl)) printf("-- err --");
            }
            fflush(stdout);

            printf("\nMy exchanges:");
            for (int r = 0; r < gt.nb_count; ++r) {
                int rank = gt.nb_ranks[r];
                printf("\n-- rank %d --", rank);

                printf("\n   snd");
                if (gt.snd_ri[r] < gt.snd_ri[r+1]) printf(" (%d .. %d)", gt.snd_ri[r], gt.snd_ri[r+1]-1);
                printf(":");
                for (int i = gt.snd_ri[r]; i < gt.snd_ri[r+1]; ++i)
                    printf(" %d", sp.l2g[gt.snd_xj[i]]);

                printf("\n   rcv");
                if (gt.rcv_ri[r] < gt.rcv_ri[r+1]) printf(" (%d .. %d)", gt.rcv_ri[r], gt.rcv_ri[r+1]-1);
                printf(":");
                for (int i = gt.rcv_ri[r]; i < gt.rcv_ri[r+1]; ++i)
                    printf(" %d", sp.l2g[i]);
            }

            printf("\n");
            fflush(stdout);
        }
        MPI_Barrier(gt.comm);
    }
}


int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    int my_rank, num_procs;
    MPI_Comm_rank(MCW, &my_rank);
    MPI_Comm_size(MCW, &num_procs);

    int Nx, Ny, Nz, k1, k2, Px, Py;
    double tol = 1e-6;

    bool valid_params = true;

    if (argc < 8) {
        valid_params = false;
        if (my_rank == 0)
            print_help();
    } else {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        Nz = atoi(argv[3]);
        k1 = atoi(argv[4]);
        k2 = atoi(argv[5]);
        Px = atoi(argv[6]);
        Py = atoi(argv[7]);
        if (argc > 8) tol = atof(argv[8]);
    }

    if (valid_params) {
        if (my_rank == 0) {
            if (Nx < 1 || Ny < 1 || Nz < 1) printf("Nx, Ny, Nz must be int > 0\n");
            if (k1 < 0 || k2 < 0) printf("k1, k2 must be int >= 0\n");
            if (Px < 1 || Py < 1) printf("Px, Py must be int > 0\n");
            if (Px > Nx || Py > Ny) printf("Px, Py must be <= Nx, Ny\n");
            if (Px * Py != num_procs) printf("Px * Py must be equal to number of processes\n");
            if (tol <= 0) printf("tol must be float > 0\n");
        }

        if (Nx < 1 || Ny < 1 || Nz < 1  || k1 < 0 || k2 < 0 || tol <= 0 ||
            Px < 1 || Py < 1 || Px > Nx || Px > Ny || Px*Py != num_procs ) {
            valid_params = false;
        }
    }

    if (!valid_params) {
        MPI_Finalize();
        return 0;
    }

    if (my_rank == 0) {
        printf("( %d x %d x %d ), ( %d / %d )\n", Nx, Ny, Nz, k1, k2);
        printf("( %d x %d ), tol = %g\n\n", Px, Py, tol);
        fflush(stdout);
    }

    Space sp;
    Gate gt;
    Matrix A;
    vector<double> b;

    {
        double t;

        if (my_rank == 0) { printf("Initiating portrait: "); fflush(stdout); }
        t = MPI_Wtime();

        Grid gr(Nx, Ny, Nz, Px, Py, my_rank);
        Topology tp(Nx, Ny, Nz, k1, k2);
        Portrait pt;
        tp.init(gr);
        init_task(gr, tp, sp, gt, pt);

        t = MPI_Wtime() - t;
        if (my_rank == 0) { printf("%f s\n", t); printf("Initiating exchanges: "); fflush(stdout); }
        t = MPI_Wtime();

        gt.init(sp, MCW);

        t = MPI_Wtime() - t;
        if (my_rank == 0) { printf("%f s\n", t); printf("Filling matrix: "); fflush(stdout); }
        t = MPI_Wtime();

        A.fill(sp, std::move(pt));

        t = MPI_Wtime() - t;
        if (my_rank == 0) { printf("%f s\n", t); printf("Filling vector: "); fflush(stdout); }
        t = MPI_Wtime();

        b = fill_vector(sp);

        t = MPI_Wtime() - t;
        if (my_rank == 0) { printf("%f s\n", t); fflush(stdout); }
    }

    Stats st;
    vector<double> x = shared::solve_cg(sp, gt, A, b, tol, true, &st);
    double x_l2 = sqrt(shared::dot(sp, gt, x, x));
    double x_min = shared::min(sp, gt, x);
    double x_max = shared::max(sp, gt, x);

    {
        const char *names[] = {"dot  ", "axpby", "SpMV ", "prod "};

        if (my_rank == 0)
            printf("Time per op |     min     |     avg     |     max     |\n");

        for (int i = 0; i < 4; ++i) {
            double min = gt.allreduce(st.times[i], MPI_MIN) / st.calls[i];
            double max = gt.allreduce(st.times[i], MPI_MAX) / st.calls[i];
            double avg = gt.allreduce(st.times[i], MPI_SUM) / num_procs / st.calls[i];

            if (my_rank == 0)
                printf("    %s   | %11.9f | %11.9f | %11.9f |\n", names[i], min, avg, max);
        }

        if (my_rank == 0)
            fflush(stdout);
    }

    if (my_rank == 0) {
        printf("\nControl values:\n");
        printf("  L2 = %f\n", x_l2);
        printf(" min = %f\n", x_min);
        printf(" max = %f\n", x_max);
        printf("\n");
        fflush(stdout);
    }

    if ( (Nx+Ny)*Nz <= 200 && num_procs <= 12)
        exchg_report(sp, gt);

    MPI_Finalize();
}
