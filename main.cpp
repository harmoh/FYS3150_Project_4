#include <random>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

// Declare functions
void initializeLattice(int nSpins, mat &spinMatrix, double &energy, double &magneticMoment);

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
    double energy = 0.0;
    double magneticMoment = 0.0;
    double temperature = 1.0;

    initializeLattice(nSpins, spinMatrix, energy, magneticMoment);

    cout << "Energy: " << energy << endl;
    cout << "Magnetic moment: " << magneticMoment << endl;

    vec energyDifference = zeros<mat>(17);
    for(int dE = -8; dE <= 8; dE += 4)
    {
        energyDifference(dE + 8) = exp(-dE/temperature);
    }

    //cout << "Energy difference:\n" << energyDifference << endl;

    // Monte Carlo cycles
    int mcCycles = 0;
    bool final = false;
    //for(int i = 0; i < mcCycles; i++)
    while(!final)
    {
        // Sweep over the lattice
        for(int y = 0; y < nSpins; y++)
        {
            for(int x = 0; x < nSpins; x++)
            {
                int ix = (int) (rand(gen) * (double) nSpins);
                int iy = (int) (rand(gen) * (double) nSpins);
                int deltaE = 2 * spinMatrix(ix,iy) *
                        (spinMatrix(pbc(ix,nSpins,-1),iy) + spinMatrix(ix,pbc(iy,nSpins,-1)) +
                         spinMatrix(pbc(ix,nSpins,1),iy) + spinMatrix(ix,pbc(iy,nSpins,1)));
                if(rand(gen) <= energyDifference(deltaE + 8))
                {
                    spinMatrix(ix,iy) *= -1.0;
                    magneticMoment += (double) 2 * spinMatrix(ix,iy);
                    energy += (double) deltaE;
                    cout << "DeltaE: " << deltaE << "\t Energy: " << energy << endl;
                }
            }
        }
        if(magneticMoment == -4) final = true;
        mcCycles++;
    }

    cout << "Energy: " << energy << endl;
    cout << "Magnetic moment: " << magneticMoment << endl;
    cout << "MC cycles: " << mcCycles << endl;

    return 0;
}

void initializeLattice(int nSpins, mat &spinMatrix, double &energy, double &magneticMoment)
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
