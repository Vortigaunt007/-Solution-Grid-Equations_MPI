#include "CommandLineArg.h"


CommandLineArg::CommandLineArg(int Nx, int Ny, int Nz, int k1, int k2, int p_x, int p_y, int Nt): nx(Nx), ny(Ny), nz(Nz), k1(k1), k2(k2), px(p_x), py(p_y), nt(Nt)
{

}

CommandLineArg::CommandLineArg(std::string filename)
{
    std::ifstream fin(filename);

    if (!fin)
        std::cout << "Grid: unable to open file for reading data" << std::endl;

    fin >> nx >> ny >> nz >> k1 >> k2 >> px >> py >> nt;

    fin.close();
}
