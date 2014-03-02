// Test program for checking the installation of libdiffpy shared library.
// Compile and run this code using:
//
//      c++ testlib.cpp -ldiffpy
//      ./a.out


#include <iostream>
#include <diffpy/srreal/ScatteringFactorTable.hpp>

int main(int argc, char* argv[])
{
    using namespace diffpy::srreal;
    ScatteringFactorTablePtr ntbl;
    ntbl = ScatteringFactorTable::createByType("neutron");
    ntbl->lookup("Pb");
    std::cout << "Installation of libdiffpy shared library works!\n";
    return 0;
}
