#include "MPI_Exchange.h"


MPI_Exchange::MPI_Exchange()
{

}

void MPI_Exchange::SendAndRecieve(MatrixCSR &M)
{
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    std::vector<int> Neighbours;

    const int sizePart = (int)M.Part.size();
    for (int i = 0; i < sizePart; i++) {
        if (world_rank != M.Part[i])
            Neighbours.push_back(M.Part[i]);
    }
    std::sort(Neighbours.begin(), Neighbours.end());
    Neighbours.erase(unique(Neighbours.begin(), Neighbours.end()), Neighbours.end() );

    int N = M.AI.size() - 1;
    for(int i = 0; i < N; i++) { // string number
        for(int j_ = M.AI[i]; j_ < M.AI[i+1]; j_++) {
            int j = M.AJ[j_]; // column number
            for (int k = 0; k < Neighbours.size(); k++)
                if (M.Part[j] == Neighbours[k]) { // is a neighbor
                  /*  if (world_rank == 0) {
                        std::cout << "Neigh = " << Neighbours[k] << " i = " << M.L2G[i] << " j = " << M.L2G[j] << "; ";
                    }*/
                    if (Recv.find(M.Part[j]) == Recv.end()) { // not found
                        std::vector<int> recv_vec;
                        recv_vec.push_back(M.L2G[j]);
                        Recv.insert(std::pair<int, std::vector<int>>(M.Part[j], recv_vec));
                    } else { // found
                        Recv.find(M.Part[j])->second.push_back(M.L2G[j]);
                    }
                    if (Send.find(M.Part[j]) == Send.end()) { // not found
                        std::vector<int> send_vec;
                        send_vec.push_back(M.L2G[i]);
                        Send.insert(std::pair<int, std::vector<int>>(M.Part[j], send_vec));
                    } else { // found
                        Send.find(M.Part[j])->second.push_back(M.L2G[i]);
                    }
                }
        }
    }

    sortSendAndRecieve(Send);
    sortSendAndRecieve(Recv);

//    SendAndRecieveLocal_to_Global(Send, M);
//    SendAndRecieveLocal_to_Global(Recv, M);

/*
    if (world_rank == 0) {
        for(std::map<int, std::vector<int> >::const_iterator it = Send.begin(); it != Send.end(); ++it) {
            std::cout << it->first << " -> ";
            for(int i = 0; i < it->second.size(); i++) {
                std::cout << it->second[i] << " ";
            }
            std::cout << std::endl;
        }
    }
    if (world_rank == 1) {
        for(std::map<int, std::vector<int> >::const_iterator it = Recv.begin(); it != Recv.end(); ++it) {
            std::cout << it->first << " -> " << world_rank << " = ";
            if (it->first == 0)
                for(int i = 0; i < it->second.size(); i++) {
                    std::cout << it->second[i] << " ";
                }
            std::cout << std::endl;
        }
    }

    if (world_rank == 2) {
        for(std::map<int, std::vector<int> >::const_iterator it = Recv.begin(); it != Recv.end(); ++it) {
            std::cout << it->first << " -> " << world_rank << " = ";
            if (it->first == 0)
                for(int i = 0; i < it->second.size(); i++) {
                    std::cout << it->second[i] << " ";
                }
            std::cout << std::endl;
        }
    }
    if (world_rank == 3) {
        for(std::map<int, std::vector<int> >::const_iterator it = Recv.begin(); it != Recv.end(); ++it) {
            std::cout << it->first << " -> " << world_rank << " = ";
            if (it->first == 0)
                for(int i = 0; i < it->second.size(); i++) {
                    std::cout << it->second[i] << " ";
                }
            std::cout << std::endl;
        }
    }
*/
    return;
}

void MPI_Exchange::sortSendAndRecieve(std::map<int, std::vector<int> > &buf)
{
    for(std::map<int, std::vector<int> >::iterator it = buf.begin(); it != buf.end(); ++it) {
        std::sort(it->second.begin(), it->second.end());
        it->second.erase(unique(it->second.begin(),it->second.end()), it->second.end() );
       // std::unique(it->second.begin(), it->second.end());
    }
}

void MPI_Exchange::SendAndRecieveLocal_to_Global(std::map<int, std::vector<int> > &buf, MatrixCSR &M)
{
    for(std::map<int, std::vector<int> >::iterator it = buf.begin(); it != buf.end(); ++it) {
        for(int i = 0; i < it->second.size(); i++) {
            it->second[i] = M.L2G[it->second[i]];
        }
    }
}
