#include <random>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;

ofstream ofile;

// Declare functions
void initializeLattice(int nSpins, mat &spinMatrix, double &energy, double &magneticMoment);
void metropolis(int nSpins, int mcCycles, double temp, vec &expectationValues);
void writeToFile(int nSpins, int mcCycles, double temp, vec expectationValues);

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

    vec expectationValues = zeros<mat>(6);

    // Monte Carlo cycles
    int mcCycles = atoi(argv[1]);

    // Declare new file name and add lattice size to file name
    string fileout = "Lattice";
    fileout.append(to_string(nSpins) + ".txt");
    ofile.open(fileout);

    clock_t time_initial = clock();

    // Start Metropolis algorithm with Monte Carlo sampling
    metropolis(nSpins, mcCycles, temp, expectationValues);
    writeToFile(nSpins, mcCycles, temp, expectationValues);

    ofile.close();

    clock_t time_final = clock();
    double elapsed_time = (time_final - time_initial) / (double) CLOCKS_PER_SEC;

    cout << "\nTime: " << elapsed_time << " seconds." << endl;

    return 0;
}

// Initalize lattice with spins, magnetic moment and energy
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

// Perform the Metropolis algorithm
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
        //cout << "exp(-dE/temperature) = " << exp(-dE/temp) << "\tfor dE = " << dE << endl;
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
        expectationValues(2) += fabs(energy);
        expectationValues(3) += magneticMoment;
        expectationValues(4) += magneticMoment * magneticMoment;
        expectationValues(5) += fabs(magneticMoment);

        //if(magneticMoment == -4) final = true;
        //mcCycles++;
    }

    cout << "Energy: " << energy << endl;
    cout << "Magnetic moment: " << magneticMoment << endl;
    cout << "MC cycles: " << mcCycles << endl;


}

void writeToFile(int nSpins, int mcCycles, double temp, vec expectationValues)
{
    // Normalization of the values
    double norm = 1.0/mcCycles;
    double normSpins = 1.0 / nSpins / nSpins;

    // Numerical values
    double expectVal_E = expectationValues(0)*norm;// / nSpins / nSpins;
    double expectVal_E2 = expectationValues(1)*norm;
    double expectVal_Eabs = expectationValues(2)*norm;
    double expectVal_M = expectationValues(3)*norm;// / nSpins / nSpins;
    double expectVal_M2 = expectationValues(4)*norm;// / nSpins / nSpins;
    double expectVal_Mabs = expectationValues(5)*norm;

    double expectVal_Cv = (expectVal_E2 - expectVal_Eabs * expectVal_Eabs) * normSpins / (temp * temp);
    double expectVal_X = (expectVal_M2 - expectVal_Mabs * expectVal_Mabs) * normSpins / temp;

    cout << "\nExpectation values, numerical: " << endl;
    //cout << "Energy: " << expectVal_E << endl;
    //cout << "Energy^2: " << expectationValues_E2 << endl;
    //cout << "|Energy|: " << expectVal_Eabs << endl;
    //cout << "Magnetic moment: " << expectVal_M << endl;
    //cout << "Magnetic moment^2: " << expectationValues_M2 << endl;
    //cout << "|Magnetic moment|: " << expectVal_Mabs << endl;

    cout << "Heat capacity: " << expectVal_Cv << endl;
    cout << "Susceptibility: " << expectVal_X << endl;

    // Analytical values
    double J = 1.0;
    double beta = 1.0;
    double Z = 4*cosh(8*J*beta) + 12;
    double expectValAnalytical_E = 32*J*sinh(8*J*beta) / Z;
    double expectValAnalytical_M2 = (32*exp(8*J*beta) + 32) / Z;
    double expectValAnalytical_Mabs = (8*exp(8*J*beta) + 16) / Z;
    double expectValAnalytical_Cv = ((256*J*cosh(8*J*beta)) / Z - expectValAnalytical_E *
                                     expectValAnalytical_E) * normSpins / (temp * temp);
    double expectValAnalytical_X = (expectValAnalytical_M2 - expectValAnalytical_Mabs *
                                    expectValAnalytical_Mabs) * normSpins / temp;
    //double Achi = 8*(exp(8.0/temp) + cosh(8.0/temp) + 3.0/2.0)/(temp*pow((cosh(8.0/temp)+3), 2));

    cout << "\nExpectation values, analytical:" << endl;
    //cout << "Energy: " << expectValAnalytical_E << endl;
    cout << "Heat capacity: " << expectValAnalytical_Cv << endl;
    cout << "Susceptibility: " << expectValAnalytical_X << endl;
    //cout << "Analytical Chi: " << Achi << endl;

    // Error between numerical and analytical
    cout << "Heat capacity (error): " << fabs(expectValAnalytical_Cv - expectVal_Cv) << endl;
    cout << "Susceptibility (error): " << fabs(expectValAnalytical_X - expectVal_X) << endl;

    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << temp;
    ofile << setw(15) << setprecision(8) << expectVal_E * normSpins;
    ofile << setw(15) << setprecision(8) << expectVal_Cv;
    ofile << setw(15) << setprecision(8) << expectVal_M;
    ofile << setw(15) << setprecision(8) << expectVal_X;

    /*
    ofile << setw(15) << setprecision(8) << E_ExpectationValues/NSpins/NSpins;
    ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
    ofile << setw(15) << setprecision(8) << M_ExpectationValues/NSpins/NSpins;
    ofile << setw(15) << setprecision(8) << Mvariance/temperature;
    ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues/NSpins/NSpins << endl;
    */
}
