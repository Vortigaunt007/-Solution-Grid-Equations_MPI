#ifndef BASICOPERATION_H
#define BASICOPERATION_H

#include "CommandLineArg.h"
#include "Vector.h"
#include "MatrixCSR.h"
#include "MPI_Exchange.h"

class BasicOperation
{
    int nt = 1;
public:
    long long  countSummAlpha = 0, countMult = 0, countDot = 0, countSolver = 0;
    long long  countAISummAlpha = 0, countAIMult = 0, countAIDot = 0;
    long long  countAISolver = 0;
    double timeSummAlpha = 0, timeMult = 0, timeDot = 0, timeUpdate = 0;

    BasicOperation();

    bool checkSizes(Vector &v1, Vector &v2);

    void resetCounters();

    void update(MatrixCSR &M, Vector &v, MPI_Exchange &K);

    double dotProduct(Vector &v1, Vector &v2, int size);
    void summAlpha(double alpha, Vector &v1_in, double beta, Vector &v2_in, Vector &out, int size);
    void mult_Matrix_Vector(MatrixCSR &, Vector &, Vector &, MPI_Exchange K);

    void TimeForBasicOperation(MatrixCSR &M, MPI_Exchange &K);
    double testDot(MatrixCSR &M);
    void testSummAlpha(MatrixCSR &M);
    void testMult(MatrixCSR &M, MPI_Exchange &K);
};

#endif // BASICOPERATION_H
