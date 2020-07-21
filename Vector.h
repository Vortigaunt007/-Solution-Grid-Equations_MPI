#ifndef VECTOR_H
#define VECTOR_H

#include "Grid.h"
#include "MatrixCSR.h"

class Vector
{
    std::vector <double> v;
    int N;

    void setSize();
public:
    Vector();
    Vector(int);

    int getSize() const;

    void Fill_Vector(MatrixCSR &M);

    double& operator[](int i); // v[] = a
    double operator[](int) const; // a = v[]
    Vector &operator =(const Vector &v);

    void saveToFile(std::string filename) const;
    void printVector() const;
};

double InitVectorElement(int i, MatrixCSR &M);

#endif // VECTOR_H
