#ifndef MPI_EXCHANGE_H
#define MPI_EXCHANGE_H

#include "CommandLineArg.h"
#include "MatrixCSR.h"
#include <map>
#include <utility>

class MPI_Exchange {
public:
    std::map<int, std::vector<int> > Send;
    std::map<int, std::vector<int> > Recv;

    void SendAndRecieve(MatrixCSR &M);
    void sortSendAndRecieve(std::map<int, std::vector<int> > &buf);
    void SendAndRecieveLocal_to_Global(std::map<int, std::vector<int> > &buf, MatrixCSR &M);

    MPI_Exchange();
};

#endif // MPI_EXCHANGE_H
