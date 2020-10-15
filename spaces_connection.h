#ifndef CONNECTION_H
#define CONNECTION_H

#include "global_entities.h"

#include <mpi.h>

#include <vector>
#include <set>
#include <map>

using std::vector;
using std::set;
using std::map;

struct Space {
    vector<int> l2g;
    map<int, int> g2l;
    int N_inner, N_interface, N_halo;
    int N_own, N_local;

    Space() : N_inner(0), N_interface(0), N_halo(0), N_own(0), N_local(0) { }
};


struct Portrait {
    vector<int> AI, AJ;
};


struct Gate {
    int my_rank;
    int num_procs;
    MPI_Comm comm;

    int nb_count;
    vector<int> nb_ranks;

    vector<int> rcv_ri;

    vector<int> snd_ri;
    vector<int> snd_xj;

    int tag;
    vector<double> buffer;
    vector<MPI_Request> req;
    vector<MPI_Status>  sts;

    Gate() : my_rank(0), num_procs(1), comm(MPI_COMM_WORLD), nb_count(0), tag(0) { }

    void waitall();

    void init(const Space& sp, MPI_Comm comm);

    void exchange(vector<double>& x);

    double allreduce(double r, MPI_Op op);
};


void init_task(
    const Grid& gr, const Topology& tp,  /* input  */
    Space& sp, Gate& gt, Portrait& pt   /* output */
);

#endif // CONNECTION_H
