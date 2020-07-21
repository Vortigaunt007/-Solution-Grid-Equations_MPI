#include "BasicOperation.h"

BasicOperation::BasicOperation()
{

}

bool BasicOperation::checkSizes(Vector &v1, Vector &v2)
{
    if(v1.getSize() != v2.getSize()) {
        std::cout << "Sizes aren't equale!" << std::endl;
        return false;
    } else
        return true;
}

void BasicOperation::resetCounters()
{
    countDot = 0;
    countMult = 0;
    countSolver = 0;
    countSummAlpha = 0;

    timeDot = 0;
    timeMult = 0;
    timeSummAlpha = 0;
    timeUpdate = 0;
}

void BasicOperation::update(MatrixCSR &M, Vector &v, MPI_Exchange &K)
{
    int world_rank;
    int sizeRequest;

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    sizeRequest = K.Send.size() + K.Recv.size();

    MPI_Request request[sizeRequest];
    MPI_Status status[sizeRequest];

    int countRequest = 0;
    double *m_free[K.Send.size()];
    double t = MPI_Wtime();

    for(std::map<int, std::vector<int> >::const_iterator it = K.Send.begin(); it != K.Send.end(); ++it) {
        std::vector<int> vec = it->second;
        int sizeVec = vec.size();
        double *vecSend = (double*)malloc(sizeof(double) * sizeVec);
        m_free[countRequest] = vecSend;
        for (int i = 0; i < sizeVec; i++) {
            int index = M.G2L[vec[i]];
            vecSend[i] = v[index];
        }
        int destRank = it->first;
      //  std::cout << destRank <<std::endl;
        MPI_Isend(vecSend, sizeVec, MPI_DOUBLE, destRank, 0, MPI_COMM_WORLD, &request[countRequest]);
        countRequest++;
    }
    int j = 0;

    double *m_recv = (double*)malloc(sizeof(double) * (M.size_Halo));
    for(std::map<int, std::vector<int> >::const_iterator it = K.Recv.begin(); it != K.Recv.end(); ++it) {
        std::vector<int> vec = it->second;
        int sizeVec = vec.size();
        int sourceRank = it->first;
        MPI_Irecv(&(m_recv[j]), sizeVec, MPI_DOUBLE, sourceRank, 0, MPI_COMM_WORLD, &request[countRequest]);
        for (int i = 0; i < sizeVec; i++) {
            int index = M.G2L[vec[i]];
            v[index] = m_recv[i];
        }

        j += sizeVec;
        countRequest++;

    }

    MPI_Waitall(countRequest, request, status);

    j = 0;
    for(std::map<int, std::vector<int> >::const_iterator it = K.Recv.begin(); it != K.Recv.end(); ++it) {
        std::vector<int> vec = it->second;
        int sizeVec = vec.size();
        for (int i = 0; i < sizeVec; i++) {
            int index = M.G2L[vec[i]];
            v[index] = m_recv[i + j];
        }
        j += sizeVec;

   }

    for (int i = 0; i < K.Send.size(); i++)
        free(m_free[i]);
    free(m_recv);

    timeUpdate += MPI_Wtime() - t;
}

void BasicOperation::mult_Matrix_Vector(MatrixCSR &M, Vector &in, Vector &out, MPI_Exchange K)
{
/*
    const int n = (int)M.AI.size() - 1;
    if(n != (in.getSize() - M.size_Halo)) {
        std::cout << "Function mult_Matrix_Vector: sizes aren't equale: " <<M.size_Own << " " << n << " " << in.getSize() << std::endl;
        return;
    }
*/
    long long count = 0;
    double t = MPI_Wtime();

    int myranksize, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &myranksize);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    for(int i = 0; i < M.size_Own; i++) { // string number
        double res = 0.0;
        for(int j_ = M.AI[i]; j_ < M.AI[i+1]; j_++) {
            int j = M.AJ[j_]; // column number
            double a_ij = M.A[j_];
            count += 2;
            res += a_ij * in[j];
        }
        out[i] = res;
    }

    timeMult += MPI_Wtime() - t;
    countAIMult += 5 * (count / 2);
    countAISolver += 5 * (count / 2);
    countMult += count;
    countSolver += count;

    return;
}

void BasicOperation::TimeForBasicOperation(MatrixCSR &M, MPI_Exchange &K)
{
    resetCounters();
    testDot(M);
    testSummAlpha(M);
    testMult(M, K);

    int myranksize, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &myranksize);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    double timeSP_VMWithUpdate = timeMult + timeUpdate;
    double maxDot;
    double maxAxpby;
    double maxSpUpdate;
    double maxSp_ELLPACK;
    MPI_Reduce(&timeSummAlpha, &maxAxpby, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeDot, &maxDot, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeSP_VMWithUpdate, &maxSpUpdate, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeMult, &maxSp_ELLPACK, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    double gAxpy = countSummAlpha/(timeSummAlpha * 1E9);
    double gDot = countDot/(timeDot * 1E9);
    double gSP_MV = countMult/(timeMult * 1E9);
    double gSP_MV_Update = countMult/((timeSP_VMWithUpdate) * 1E9);
    MPI_Allreduce(&gAxpy, &gAxpy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&gDot, &gDot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&gSP_MV, &gSP_MV, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&gSP_MV_Update, &gSP_MV_Update, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(world_rank == 0) {
        std::cout << "_______________________" << std::endl;
        std::cout << "Time AXPBY operation: " << maxAxpby << std::endl;
        std::cout << "Time SpMV_ELLPACK operation: " << maxSp_ELLPACK << std::endl;
        std::cout << "Time SpMV_with_update: operation:" << maxSpUpdate << std::endl;
        std::cout << "Time DOT: operation:" << maxDot << std::endl;
        std::cout << "SpMV_ELLPACK GFLOPS: " << gSP_MV << std::endl;
        std::cout << "SpMV_ELLPACK_UPDATE GFLOPS: " << gSP_MV_Update << std::endl;
        std::cout << "AXPBY GFLOPS: " << gAxpy << std::endl;
        std::cout << "DOT usually GFLOPS: " << gDot << std::endl;
        std::cout << "_______________________" << std::endl;
    }
}


double BasicOperation::testDot(MatrixCSR &M)
{
    Vector a(M.size_Own);
    Vector b(M.size_Own);
    a.Fill_Vector(M);
    b.Fill_Vector(M);

    double res;
    for (int i = 0; i < 150; i++)
        res = dotProduct(a, b, M.size_Own);
    return res;
}

void BasicOperation::testSummAlpha(MatrixCSR &M)
{
    Vector a(M.size_Own);
    Vector b(M.size_Own);
    a.Fill_Vector(M);
    b.Fill_Vector(M);

    Vector c(M.size_Own);

    for (int i = 0; i < 150; i++)
        summAlpha(2.0, a, 4.0, b, c, M.size_Own);

    return;
}

void BasicOperation::testMult(MatrixCSR &M, MPI_Exchange &K)
{
    Vector a(M.size_Own + M.size_Halo);
    Vector b(M.size_Own + M.size_Halo);
    a.Fill_Vector(M);

    for (int i = 0; i < 150; i++) {
        update(M, a, K);
        mult_Matrix_Vector(M, a, b, K);
    }
    return;
}

void BasicOperation::summAlpha(double alpha, Vector &v1_in, double beta, Vector &v2_in, Vector &out, int size)
{
   // if(!checkSizes(v1_in, v2_in))
   //     return;

    double t = MPI_Wtime();

    for(int i = 0; i < size; i++) {
        out[i] = alpha * v1_in[i] + beta * v2_in[i];
    }

    timeSummAlpha += MPI_Wtime() - t;
    countAISummAlpha += 5 * size;
    countAISolver += 5 * size;
    countSolver += 3 * size;
    countSummAlpha += 3 * size;

    return;
}

double BasicOperation::dotProduct(Vector &v1, Vector &v2, int size)
{
    if(!checkSizes(v1, v2))
        return 0.0;

    double result = 0.0;
    double t = MPI_Wtime();

    for(int i = 0; i < size; i++) {
        result += v1[i] * v2[i];
    }

    double resDotProductAll;
    MPI_Allreduce(&result, &resDotProductAll, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    timeDot += MPI_Wtime() - t;    
    countAIDot += 4 * size;
    countAISolver += 4 * size;
    countSolver += 2 * size;
    countDot += 2 * size;

    return resDotProductAll;
}
