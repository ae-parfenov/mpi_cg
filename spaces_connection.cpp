#include "spaces_connection.h"

void Gate::waitall() {
    if (nb_count == 0) return;
    MPI_Waitall(nb_count*2, req.data(), sts.data());
}

void Gate::init(const Space& sp, MPI_Comm comm) {

    this->comm = comm;
    tag = 0;
    req.resize(2*nb_count);
    sts.resize(2*nb_count);

    vector<int> snd_buf(nb_count);
    vector<int> rcv_buf(nb_count);
    int k = 0;
    for (int i = 0; i < nb_count; ++i) {
        snd_buf[i] = rcv_ri[i+1] - rcv_ri[i];
        MPI_Isend(&snd_buf[i], 1, MPI_INT, nb_ranks[i], tag, comm, &req[k++]);
        MPI_Irecv(&rcv_buf[i], 1, MPI_INT, nb_ranks[i], tag, comm, &req[k++]);
    }

    waitall();
    int snd_count = 0;
    snd_ri.push_back(0);
    for (int i = 0; i < nb_count; ++i) {
        snd_count += rcv_buf[i];
        snd_ri.push_back(snd_count);
    }

    snd_xj.resize(snd_count);
    buffer.resize(snd_count);

    ++tag;
    k = 0;
    for (int i = 0; i < nb_count; ++i) {
        MPI_Isend(&sp.l2g[rcv_ri[i]], rcv_ri[i+1]-rcv_ri[i], MPI_INT, nb_ranks[i], tag, comm, &req[k++]);
        MPI_Irecv(&snd_xj[snd_ri[i]], snd_ri[i+1]-snd_ri[i], MPI_INT, nb_ranks[i], tag, comm, &req[k++]);
    }

    waitall();
    for (int i = 0; i < snd_count; ++i)
        snd_xj[i] = sp.g2l.at(snd_xj[i]);
}

void Gate::exchange(vector<double>& x) {
    ++tag;
    if (nb_count == 0) return;
    int k = 0;

    for (int i = 0; i < nb_count; ++i) {
        for (int j = snd_ri[i]; j < snd_ri[i+1]; ++j)
            buffer[j] = x[snd_xj[j]];
        MPI_Isend(&buffer[snd_ri[i]], snd_ri[i+1]-snd_ri[i], MPI_DOUBLE, nb_ranks[i], tag, comm, &req[k++]);
        MPI_Irecv(&x[rcv_ri[i]], rcv_ri[i+1]-rcv_ri[i], MPI_DOUBLE, nb_ranks[i], tag, comm, &req[k++]);
    }
}

double Gate::allreduce(double r, MPI_Op op) {
    double z = r;
    if (num_procs > 1)
        MPI_Allreduce(&r, &z, 1, MPI_DOUBLE, op, comm);
    return z;
}


void init_task(
    const Grid& gr, const Topology& tp,  /* input  */
    Space& sp, Gate& gt, Portrait& pt   /* output */
) {

    int N_own = tp.AI.size()-1;
    int my_id = gr.id;
    vector<int> nbs = gr.neighboring_processes();
    int nb_count = nbs.size();

    vector<int> inner_nodes;
    vector<int> interface_nodes;

    map<int, set<int> > halo_nodes;
    for (int r = 0; r < nb_count; ++r)
        halo_nodes[nbs[r]] = set<int>();

    for (int I = 0; I < N_own; ++I) {
        bool is_interface = false;
        for (int J = tp.AI[I]; J < tp.AI[I+1]; ++J) {
            int gl_index = tp.AJ[J];
            int rank = gr.part(gl_index);
            if (rank != my_id) {
                is_interface = true;
                if (halo_nodes[rank].count(gl_index) == 0)
                    halo_nodes[rank].insert(gl_index);
            }
        }
        if (is_interface)
            interface_nodes.push_back(I);
        else
            inner_nodes.push_back(I);
    }

    int N_inner = sp.N_inner = inner_nodes.size();
    sp.N_interface = interface_nodes.size();
    sp.N_own = N_own;

    int I = 0;
    for (int i = 0; i < N_own; ++i) {
        int index = i < N_inner ? inner_nodes[i] : interface_nodes[i-N_inner];
        int gl_index = tp.gl_ids[index];
        sp.l2g.push_back(gl_index);
        sp.g2l[gl_index] = I++;
    }

    gt.my_rank = gr.id;
    gt.num_procs = gr.Px * gr.Py;
    gt.nb_count = nb_count;

    for (int r = 0; r < nb_count; ++r) {
        int rank = nbs[r];
        gt.rcv_ri.push_back(I);
        for (
                set<int>::const_iterator It = halo_nodes[rank].cbegin();
                It != halo_nodes[rank].cend();
                It++
            ) {
            int gl_index = *It;
            sp.l2g.push_back(gl_index);
            sp.g2l[gl_index] = I++;
        }
    }

    gt.rcv_ri.push_back(I);
    gt.nb_ranks = std::move(nbs);

    sp.N_local = sp.l2g.size();
    sp.N_halo = sp.N_local - sp.N_own;

    pt.AI.push_back(0);
    for (int i = 0; i < N_own; ++i) {
        int index = i < N_inner ? inner_nodes[i] : interface_nodes[i-N_inner];
        set<int> nodes;
        for (int j = tp.AI[index]; j < tp.AI[index+1]; ++j)
            nodes.insert(sp.g2l[tp.AJ[j]]);
        for (set<int>::const_iterator It = nodes.cbegin(); It != nodes.cend(); It++)
            pt.AJ.push_back(*It);
        pt.AI.push_back(pt.AJ.size());
    }
}
