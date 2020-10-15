#ifndef GLOBAL_H
#define GLOBAL_H

#include <vector>

using std::vector;

struct Grid {
    int Nx, Ny, Nz;
    int Px, Py;
    int id;

    int i0, i1;
    int j0, j1;
    int k0, k1;

    Grid(int Nx, int Ny, int Nz, int Px, int Py, int id);

    vector<int> neighboring_processes() const;

    int part(int gl_id) const;
};


struct Topology {
    int Nx, Ny, Nz;
    int K1, K2;

    vector<int> gl_ids;
    vector<int> AI, AJ;

    Topology(int Nx, int Ny, int Nz, int K1, int K2)
    : Nx(Nx), Ny(Ny), Nz(Nz), K1(K1), K2(K2) { };

    int global_id(int i, int j, int k) const;

    bool is_inside(int i, int j, int k) const;

    int add_adjacent(vector<int> &nbs, int i, int j, int k) const;


    void init(const Grid &gr);
};

#endif // GLOBAL_H
