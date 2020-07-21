#include "Grid.h"


int Grid::getPx() const
{
    return px;
}

int Grid::getPy() const
{
    return py;
}

Grid::Grid(): Nx(2), Ny(2), Nz(2), k1(1), k2(1)
{

}

Grid::Grid(int x, int y, int z, int k_1, int k_2, int p_x, int p_y): Nx(x), Ny(y), Nz(z), k1(k_1), k2(k_2), px(p_x), py(p_y)
{

}

Grid::Grid(std::string filename)
{
    std::ifstream fin(filename);

    if (!fin)
        std::cout << "Grid: unable to open file for reading data" << std::endl;

    fin >> Nx >> Ny >> Nz >> k1 >> k2;

    fin.close();
}

int Grid::getIndex(int i, int j, int k)
{
    return (i + j * Nx + k * Nx * Ny);
}

std::vector<int> Grid::getCoordinates(int index)
{
    std::vector<int> coord(3);

    int i, j, k;

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    k = index / (getNx() * getNy());
    j = (index / getNx()) % getNy();
    i = index % getNx();

    coord[0] = i;
    coord[1] = j;
    coord[2] = k;

    return coord;
}

cell_type Grid::TypeCell(int i, int, int k)
{
    if (k == 0 || (k % (k1 + k2) == 0 && i == 0) || (k % (k1 + k2) == 1 && i == (Nx - 1))) // border && one layer of nodes at the beginning && end of the prism layer
        return HEX;
    if (k1 > 1 && (0 < k % (k1 + k2)) && (k % (k1 + k2) < k1)) // nodes inside the hexahedron layer
        return HEX;
    else
        return PRISM;
}

std::vector<int> Grid::Neighboring_Nodes(int i, int j, int k, int begin_x, int end_x, int begin_y, int end_y)
{
    std::vector <int> neighb_N;

    // robot with one step
    int dx = 1, dy = 0, dz = 0;
    int nx = 0, ny = 0, nz = 0;

    neighb_N.push_back(getIndex(i, j, k)); // the node itself
    for(int l = 1; l < 7; l++) {
        nx = i + dx;
        ny = j + dy;
        nz = k + dz;

        if ((0 <= nx) && (nx < Nx) && (0 <= ny) && (ny < Ny) && (0 <= nz) && (nz < Nz)) {
            neighb_N.push_back(getIndex(nx, ny, nz));
        }

        if (l == 1) {
            dx = -1; dy = 0; dz = 0;
        }
        else if (l == 2) {
            dx = 0; dy = 1; dz = 0;
        }
        else if (l == 3) {
            dx = 0; dy = -1; dz = 0;
        }
        else if (l == 4) {
            dx = 0; dy = 0; dz = 1;
        }
        else if (l == 5) {
            dx = 0; dy = 0; dz = -1;
        }
    }

    if (TypeCell(i, j, k) == PRISM) {
        if ((k % (k1 + k2) == 0) && ((i-1) < Nx) && ((i-1) >= 0) && ((k-1) < Nz) && ((k-1) >= 0))  // border end
            neighb_N.push_back(getIndex(i - 1, j, k - 1));
        else if (k % (k1 + k2) == k1 && ((i+1) < Nx) && ((i+1) >= 0) && ((k+1) < Nz) && ((k+1) >= 0)) // border begin
            neighb_N.push_back(getIndex(i + 1, j, k + 1));
        else  {// internal nodes of the layer
            if ((i != Nx - 1) && (k != Nz - 1) && ((i+1) < Nx) && ((i+1) >= 0) && ((k+1) < Nz) && ((k+1) >= 0))
                neighb_N.push_back(getIndex(i + 1, j, k + 1));
            if ((i != 0) && (k != 0) && ((i-1) < Nx) && ((i-1) >= 0) && ((k-1) < Nz) && ((k-1) >= 0))
                neighb_N.push_back(getIndex(i - 1, j, k - 1));
        }
    }

    sort(neighb_N.begin(), neighb_N.end());

    return neighb_N;
}


int Grid::getK2() const
{
    return k2;
}

int Grid::getNx() const
{
    return Nx;
}

int Grid::getNy() const
{
    return Ny;
}

int Grid::getNz() const
{
    return Nz;
}


int Grid::getK1() const
{
    return k1;
}
