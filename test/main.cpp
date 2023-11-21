#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include <fenv.h>

#include "Field.hpp"
#include "Vars.hpp"
#include "Mat.hpp"


int main(int argc, char** argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    Mat<3,3> A = Mat<3,3>();

    A[0][0] = 1;
    A[0][1] = 2;
    A[0][2] = 3;
    A[1][0] = 5;
    A[1][1] = 4;
    A[1][2] = 1;
    A[2][0] = 6;
    A[2][1] = 7;
    A[2][2] = 4;

    Mat<3,3> C = A;

    C[0][2]= 5;

    Field<Mat<3,3>> matrixes = Field<Mat<3,3>>(2);

    matrixes[0] = A;
    matrixes[1] = C;

    Vars<3>b = Vars<3>({1,2,3});
    Vars<3>d = Vars<3>({3,2,1});

    Field<Vars<3>> rhs = Field<Vars<3>>(2);

    rhs[0] = b;
    rhs[1] = d;


    Field<Vars<3>> res = fdot(matrixes, rhs);

    return 0;
}