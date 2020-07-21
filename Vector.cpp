#include "Vector.h"

Vector::Vector(): N(10)
{
    setSize();
}

Vector::Vector(int str): N(str)
{
    setSize();
}

void Vector::setSize()
{
    v.resize(N);
}

int Vector::getSize() const
{
    return N;
}

void Vector::Fill_Vector(MatrixCSR &M)
{
    for(int i = 0; i < N; i++)
        v[i] = cos(double(M.L2G[i]) * double(M.L2G[i]));;
}

double &Vector::operator[](int i)
{
      return v[i];
}

double Vector::operator[](int i) const
{

    return v[i];
}

Vector &Vector::operator =(const Vector &v2)
{
    for (int i = 0; i < N; i++)
        v[i] = v2[i];
    return *this;
}


double InitVectorElement(int i, MatrixCSR &M)
{
    return cos(double(M.L2G[i]) * double(M.L2G[i]));
}



void Vector::saveToFile(std::string filename) const
{
    std::ofstream fout(filename);

    for(int i = 0; i < N; i++)
        fout << v[i] << " ";
    fout << std::endl;

    fout.close();
}

void Vector::printVector() const
{
    std::cout << "Vector v" << std::endl;

    for(int i = 0; i < N; i++)
        std::cout << v[i] << ' ';
    std::cout << std::endl;
}
