#ifndef MATRIXCSR_H
#define MATRIXCSR_H

#include "Grid.h"

class MatrixCSR
{
public:
    std::vector<int> AI;
    std::vector<int> AJ;
    std::vector<double> A;

    std::vector<int> L2G;
    std::vector<int> G2L;
    std::vector<int> Part;

    int size_Halo;
    int size_Own;

    std::vector<std::vector<int>> range_table;

    MatrixCSR();
    MatrixCSR(std::vector<int> &AI, std::vector<int> &AJ, std::vector<double> &A);

    void Matrix_Portrait(Grid);
    void Fill_Matrix(double f(int, int));

    void saveToFile(std::string filename) const;
    void printMatrixCSR_Size() const;
    void printMatrixCSR_A() const;
    
    std::vector<int> getCoordForProc(Grid D);
    int getProcForCoord(Grid D, int proc);

    int IntervalBegin(int n, int p, int proc_x);
    int IntervalEnd(int n, int p, int proc_x);
};

void inverse_Matrix(MatrixCSR &A, MatrixCSR &A_1);
double InitDiagElement(double);
double InitOffDiagElement(int i, int j);

#endif // MATRIXCSR_H
