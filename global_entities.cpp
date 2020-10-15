#include "global_entities.h"

Grid::Grid(int Nx, int Ny, int Nz, int Px, int Py, int id)
    : Nx(Nx), Ny(Ny), Nz(Nz), Px(Px), Py(Py), id(id) {
    int Pi = id % Px;
    int Pj = id / Px;

    i0 = (Nx / Px) * Pi + (Pi <= Nx % Px ? Pi : Nx % Px);
    j0 = (Ny / Py) * Pj + (Pj <= Ny % Py ? Pj : Ny % Py);

    i1 = i0 + (Nx / Px) + int(Nx % Px > Pi);
    j1 = j0 + (Ny / Py) + int(Ny % Py > Pj);

    k0 = 0;
    k1 = Nz;
}

vector<int> Grid::neighboring_processes() const {
    int Pi = id % Px;
    int Pj = id / Px;

    vector<int> nbs;
    if (Pj > 0) nbs.push_back(id - Px);
    if (Pi > 0) nbs.push_back(id - 1);
    if (Pi < Px-1) nbs.push_back(id + 1);
    if (Pj < Py-1) nbs.push_back(id + Px);

    return nbs;
}

int Grid::part(int gl_id) const {
    gl_id = gl_id % (Nx*Ny);
    int i = gl_id % Nx;
    int j = gl_id / Nx;

    int Kx = Nx / Px;
    int Dx = Nx % Px;
    int Pi = i < Dx*(Kx+1) ? i / (Kx+1) : (i - Dx*(Kx+1)) / Kx + Dx;

    int Ky = Ny / Py;
    int Dy = Ny % Py;
    int Pj = j < Dy*(Ky+1) ? j / (Ky+1) : (j - Dy*(Ky+1)) / Ky + Dy;

    return Pj*Px + Pi;
}


int Topology::global_id(int i, int j, int k) const {
    return Nx*Ny*k + Nx*j + i;
}

bool Topology::is_inside(int i, int j, int k) const {
    return (0 <= i && i < Nx && 0 <= j && j < Ny && 0 <= k && k < Nz);
}

int Topology::add_adjacent(vector<int> &nbs, int i, int j, int k) const {
    static int di[] = {0, 0, -1, 0, 1, 0, 0};
    static int dj[] = {0, -1, 0, 0, 0, 1, 0};
    static int dk[] = {-1, 0, 0, 0, 0, 0, 1};

    int size = nbs.size();
    bool prism_prev = (K2 > 0) && (k > 0) && ( (k-1) % (K1+K2) >= K1 );
    bool prism_next = (K2 > 0) && (k < Nz-1) && ( k % (K1+K2) >= K1 );

    if (prism_prev && is_inside(i-1, j, k-1))
        nbs.push_back(global_id(i-1, j, k-1));

    for (int l = 0; l < 7; ++l)
        if (is_inside(i+di[l], j+dj[l], k+dk[l]))
            nbs.push_back(global_id(i+di[l], j+dj[l], k+dk[l]));

    if (prism_next && is_inside(i+1, j, k+1))
        nbs.push_back(global_id(i+1, j, k+1));

    return nbs.size() - size;
}

void Topology::init(const Grid &gr) {
    gl_ids.resize(0);
    AI.resize(0);
    AJ.resize(0);

    AI.push_back(0);
    for (int k = gr.k0; k < gr.k1; ++k)
        for (int j = gr.j0; j < gr.j1; ++j)
            for (int i = gr.i0; i < gr.i1; ++i) {
                int id = global_id(i,j,k);
                add_adjacent(AJ, i,j,k);
                AI.push_back(AJ.size());
                gl_ids.push_back(id);
            }
}
