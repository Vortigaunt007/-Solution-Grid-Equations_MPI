#include "CommandLineArg.h"
#include "Grid.h"
#include "MatrixCSR.h"
#include "Vector.h"
#include "Solver.h"
#include "MPI_Exchange.h"

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    CommandLineArg arg("input.txt"); // Nx, Ny, Nz, k1, k2, px, py, nt

    Grid D(arg.nx, arg.ny, arg.nz, arg.k1, arg.k2, arg.px, arg.py);

    MatrixCSR M;
    M.Matrix_Portrait(D);
    Vector b(M.size_Own + M.size_Halo);

    MPI_Exchange Buf_mpi;
    Buf_mpi.SendAndRecieve(M);

    Solver S(M, b);
    double res = S.ConjugateGradient(Buf_mpi);

    BasicOperation B;
    B.TimeForBasicOperation(M, Buf_mpi);

    MPI_Finalize();
    return 0;
}
