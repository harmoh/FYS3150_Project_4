#include <random>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

// Declare functions
void initializeLattice(int nSpins, mat spinMatrix, double &energy, double &magneticMoment);

// Periodic boundary conditions
int pbc(int i, int limit, int add)
{
    return (i + limit + add) % limit;
}

// Main program
int main(int argc, char *argv[])
{
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> rand(0.0,1.0);

    int nSpins = 2;
    mat spinMatrix = zeros<mat>(nSpins,nSpins);
    double energy = 0;
    double magneticMoment = 0;

    cout << "Energy: " << energy << endl;
    cout << "Magnetic moment: " << magneticMoment << endl;

    initializeLattice(nSpins, spinMatrix, energy, magneticMoment);

    //cout << "Matrix:\n" << spinMatrix << endl;
    cout << "Energy: " << energy << endl;
    cout << "Magnetic moment: " << magneticMoment << endl;

    return 0;
}

void initializeLattice(int nSpins, mat spinMatrix, double &energy, double &magneticMoment)
{
    // Set ground state of the lattice and magnetic moment
    for(int y = 0; y < nSpins; y++)
    {
        for(int x = 0; x < nSpins; x++)
        {
            spinMatrix(x,y) = 1.0;
            magneticMoment += (double) spinMatrix(x,y);
        }
    }

    // Set up energy
    for(int y = 0; y < nSpins; y++)
    {
        for(int x = 0; x < nSpins; x++)
        {
            energy -= spinMatrix(x,y)*(spinMatrix(pbc(x,nSpins,-1),y) + spinMatrix(x,pbc(y,nSpins,-1)));
        }
    }
}
