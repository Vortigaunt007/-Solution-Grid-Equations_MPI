#include "MatrixCSR.h"


MatrixCSR::MatrixCSR()
{

}

MatrixCSR::MatrixCSR(std::vector<int> &AI_t, std::vector<int> &AJ_t, std::vector<double> &A_t)
{
    const int n = (int)A_t.size();
    for(int i = 0; i < n; i++)
        A.push_back(A_t[i]);

    const int m = AI_t.size();
    for(int i = 0; i < m; i++)
        AI.push_back(AI_t[i]);

    const int k = AJ_t.size();
    for(int i = 0; i < k; i++)
        AJ.push_back(AJ_t[i]);
}

void MatrixCSR::Matrix_Portrait(Grid D)
{
    G2L.resize(D.getNx() * D.getNy() * D.getNz());
    std::vector<int> sort_list_of_neighbors;
    int size = 0;

    int myranksize, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &myranksize);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::vector<int> proc_mpi = getCoordForProc(D);

    int begin_x, end_x, begin_y, end_y;
    begin_x = IntervalBegin(D.getNx(), D.getPx(), proc_mpi[0]);
    end_x = IntervalEnd(D.getNx(), D.getPx(), proc_mpi[0]);
    begin_y = IntervalBegin(D.getNy(), D.getPy(), proc_mpi[1]);
    end_y = IntervalEnd(D.getNy(), D.getPy(), proc_mpi[1]);

    for(int k = 0; k < D.getNz(); k++)
        for(int j = begin_y; j < end_y; j++)
            for(int i = begin_x; i < end_x; i++) {
                //std::cout <<  "hey, my rank" << world_rank << ", index = "<<  D.getIndex(i, j, k) << " "<<i<<" "<<j<<" "<< k<< std::endl;
                L2G.push_back(D.getIndex(i, j, k));
                sort_list_of_neighbors = D.Neighboring_Nodes(i, j, k, begin_x, end_x, begin_y, end_y); // list of neighboring nodes for a node with coordinates (i, j, k)
                size = AJ.size();
                AI.push_back(size);
                const int bound = sort_list_of_neighbors.size();
                for(int l = 0; l < bound; l++) {
                    int neighbor_index = sort_list_of_neighbors[l];
                    /*if (world_rank == 0) {
                        std::cout << neighbor_index << std::endl;
                    }*/
                    AJ.push_back(neighbor_index);
                   // std::cout <<  "hey, my rank" << world_rank << ", neighb = "<<  neighbor_index << std::endl;
                    A.push_back(1.0); // random value
                }
            }

    size = AJ.size();
    AI.push_back(size);

    const int mc = (int)G2L.size();
    for (int i = 0; i < mc; i++)
        G2L[i] = -1;

    const int nc = (int)L2G.size();
    for (int i = 0; i < nc; i++)
        G2L[L2G[i]] = i;

    size_Own = AI.size() - 1;

    int ncount = AI.size() - 1; // N_0
    const int count = (int)AJ.size();
    for (int i = 0; i < count; i++) {
        if (G2L[AJ[i]] < 0) {
            L2G.push_back(AJ[i]); // L2G[ncount] = AJ[i]
            G2L[AJ[i]] = ncount; // G2L[AJ[i]] = ncount
            ncount++;
        }
    }

    size_Halo = L2G.size() - size_Own;

    int val;
    for (int i = 0; i < count; i++) {
        val = AJ[i];
        AJ[i] = G2L[val];
    }

    const int size_L2G = (int)L2G.size();
    int proc_num;

    for(int i = 0; i < size_L2G; i++) {
        proc_num = getProcForCoord(D, L2G[i]);
        Part.push_back(proc_num);
    }
}

void MatrixCSR::Fill_Matrix(double f(int, int))
{
    int diag_elem_index = 0;
    double elem = 0;
    double sum_of_row_elem = 0.0; // sum of row elements except diagonal

    int myranksize, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &myranksize);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int N = AI.size() - 1;

    for(int i = 0; i < N; i++) { // string number
        sum_of_row_elem = 0.0;
        for(int j_ = AI[i]; j_ < AI[i+1]; j_++) {
            int j = L2G[AJ[j_]]; // column number
            if (L2G[i] == j) { // diagonal element
                diag_elem_index = j_; // we memorize the index in matrix A so that then we push the desired value there
                A[j_] = 1.0; // random value
            }
            else {
                elem = f(L2G[i], j);
                A[j_] = elem;
                sum_of_row_elem += abs(elem);
            }
        }
        A[diag_elem_index] = InitDiagElement(sum_of_row_elem); // diagonal element
    }
}


void inverse_Matrix(MatrixCSR &A, MatrixCSR &A_1)
{
    const int n = (int)A.AI.size();
    for(int i = 0; i < n; i++)
        A_1.AI.push_back(A.AI[i]);

    const int m = (int)A.AJ.size();
    for(int i = 0; i < m; i++)
        A_1.AJ.push_back(A.AI[i]);

    const int k = (int)A.A.size();
    for(int i = 0; i < k; i++)
        A_1.A.push_back(1.0 / A.A[i]); // diagonal matrix
}

double InitDiagElement(double a)
{
    return 1.5 * a;
}

double InitOffDiagElement(int i, int j)
{
    return cos(double(i) * double(j));
    //return sin(i + j + 1);
}

void MatrixCSR::saveToFile(std::string filename) const
{
    std::ofstream fout(filename);

    const int n = (int)A.size();
    for(int i = 0; i < n; i++)
        fout << A[i] << " ";
    fout << std::endl;

    const int m = (int)AI.size();
    for(int i = 0; i < m; i++)
        fout << AI[i] << " ";
    fout << std::endl;

    const int k = (int)AJ.size();
    for(int i = 0; i < k; i++)
        fout << AJ[i] << " ";
    fout << std::endl;

    fout.close();
}

void MatrixCSR::printMatrixCSR_A() const
{
    std::cout << "MatrixCSR A.A" << std::endl;

    const int n = (int)A.size();
    for(int i = 0; i < n; i++)
        std::cout << A[i] << " ";
    std::cout<< std::endl;
}

std::vector<int> MatrixCSR::getCoordForProc(Grid D)
{
    std::vector<int> coord(2);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    coord[0] = world_rank % D.getPx();
    coord[1] = world_rank / D.getPx();

    return coord;
}

int MatrixCSR::getProcForCoord(Grid D, int glob_num)
{
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // process grid range table
    range_table.resize(world_size);
    const int table_size = (int)range_table.size();
    for (int i = 0; i < table_size; i++) {
        range_table[i].resize(4); // [begin_x, end_x, begin_y, end_y]
    }

    for (int i = 0; i < table_size; i++) {
        range_table[i][0] = IntervalBegin(D.getNx(), D.getPx(), i % D.getPx());
        range_table[i][1] = IntervalEnd(D.getNx(), D.getPx(), i % D.getPx());
        range_table[i][2] = IntervalBegin(D.getNy(), D.getPy(), i / D.getPx());
        range_table[i][3] = IntervalEnd(D.getNy(), D.getPy(), i / D.getPx());
    }

    std::vector<int> v(3);
    v = D.getCoordinates(glob_num);

    int num_proc;
    for (int i = 0; i < table_size; i++) {
        if ((v[0] >= range_table[i][0]) && (v[0] < range_table[i][1])
                && (v[1] >= range_table[i][2]) && (v[1] < range_table[i][3])) {
            num_proc = i;
            break;
        }
    }

    return num_proc;
 }

int MatrixCSR::IntervalBegin(int n, int p, int pos_procGrid)
{
    int myranksize;
    int world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &myranksize);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    return pos_procGrid * (p > 0 ? n / p : 0) + std::min(pos_procGrid, n % p);
}

int MatrixCSR::IntervalEnd(int n, int p, int pos_procGrid)
{
    int myranksize;
    int world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &myranksize); // number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // current process number

    int begin = pos_procGrid * (n / p) + std::min(pos_procGrid, n % p);

    return begin + (p > 0 ? n / p : 0) + (pos_procGrid < (n % p));
}

void MatrixCSR::printMatrixCSR_Size() const
{
    std::cout << "A.size = " << A.size() << std::endl;
    std::cout << "AI.size = " << AI.size() << std::endl;
    std::cout << "AJ.size = " << AJ.size() << std::endl;
}
