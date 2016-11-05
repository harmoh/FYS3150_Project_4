#include <random>
#include <iostream>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;

// Declare functions
void initializeLattice(int nSpins, mat &spinMatrix, double &energy, double &magneticMoment);
void metropolis(int nSpins, int mcCycles, double temp, vec &expectationValues);

// Periodic boundary conditions
int pbc(int i, int limit, int add)
{
    return (i + limit + add) % limit;
}

// Main program
int main(int argc, char *argv[])
{
    int nSpins = 2;
    double temp = 1.0;

    //cout << "Energy: " << energy << endl;
    //cout << "Magnetic moment: " << magneticMoment << endl;

    //cout << "Energy difference:\n" << energyDifference << endl;

    vec expectationValues = zeros<mat>(5);

    // Monte Carlo cycles
    int mcCycles = atoi(argv[1]);
    clock_t time_initial = clock();
    metropolis(nSpins, mcCycles, temp, expectationValues);
    clock_t time_final = clock();
    double elapsed_time = (time_final - time_initial) / (double) CLOCKS_PER_SEC;

    cout << "Time: " << elapsed_time << " seconds." << endl;

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

void metropolis(int nSpins, int mcCycles, double temp, vec &expectationValues)
{
    // Initialize RNG, can be called by rand(gen) to get a random number between 0 and 1
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> rand(0.0,1.0);

    double energy = 0.0;
    double magneticMoment = 0.0;
    mat spinMatrix = zeros<mat>(nSpins,nSpins);
    initializeLattice(nSpins, spinMatrix, energy, magneticMoment);

    vec energyDifference = zeros<mat>(17);
    for(int dE = -8; dE <= 8; dE += 4)
    {
        energyDifference(dE + 8) = exp(-dE/temp);
        cout << "exp(-dE/temperature) = " << exp(-dE/temp) << "\tfor dE = " << dE << endl;
    }

    for(int i = 0; i < mcCycles; i++)
    //while(!final)
    {
        // Sweep over the lattice
        for(int x = 0; x < nSpins; x++)
        {
            for(int y = 0; y < nSpins; y++)
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
                    //cout << "DeltaE: " << deltaE << "\t Energy: " << energy << endl;
                }
            }
        }

        expectationValues(0) += energy;
        expectationValues(1) += energy * energy;
        expectationValues(2) += magneticMoment;
        expectationValues(3) += magneticMoment * magneticMoment;
        expectationValues(4) += fabs(magneticMoment);

        //if(magneticMoment == -4) final = true;
        //mcCycles++;
    }

    cout << "Energy: " << energy << endl;
    cout << "Magnetic moment: " << magneticMoment << endl;
    cout << "MC cycles: " << mcCycles << endl;

    double norm = 1.0/mcCycles;
    double E_ExpectationValues = expectationValues(0)*norm;
    double E2_ExpectationValues = expectationValues(1)*norm;
    double M_ExpectationValues = expectationValues(2)*norm;
    double M2_ExpectationValues = expectationValues(3)*norm;
    double Mabs_ExpectationValues = expectationValues(4)*norm;

    cout << "Expectation values:\nEnergy: " << E_ExpectationValues/nSpins/nSpins << "\nEnergy^2: " <<
            E2_ExpectationValues << "\nMagnetic moment: " << M_ExpectationValues/nSpins/nSpins <<
            "\nMagnetic moment^2: " << M2_ExpectationValues << "\n|Magnetic moment|: " <<
            Mabs_ExpectationValues << endl;
}
