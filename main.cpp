#include "mpi.h"
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
void metropolis(int nSpins, mat &spinMatrix, double &energy, double &magneticMoment, double &acceptedFlips, vec &energyDifference);
void writeToFile(int nSpins, double mcCycles, double temp, vec expectationValues, double &CvError, double &XError, double acceptedFlips);

// Periodic boundary conditions
int pbc(int i, int limit, int add)
{
    return (i + limit + add) % limit;
}

// Main program
int main(int argc, char *argv[])
{
    int nSpins, numberOfLoops, my_rank, numprocs;
    double tempInit, tempFinal, tempStep, mcCycles;

    //  MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    // Get nSpins, mcCycles, tempInit, tempFinal and tempStep from input arguments
    if (argc < 6)
    {
        cout << "Bad Usage: " << argv[0] << endl;
        cout << "Read number of spins, MC cycles, initial and final temperature, ";
        cout << "tempurate step and number of loops." << endl;
        exit(1);
    }
    else
    {
        nSpins = atoi(argv[1]);
        mcCycles = atof(argv[2]);
        tempInit = atof(argv[3]);
        tempFinal = atof(argv[4]);
        tempStep = atof(argv[5]);
        numberOfLoops = atoi(argv[6]);
    }
    if(my_rank == 0)
    {
        cout << "MC cycles: " << mcCycles << endl;
    }

    /*
     * Determine number of intervall which are used by all processes
     * myloop_begin gives the starting point on process my_rank
     * myloop_end gives the end point for summation on process my_rank
     */
    int noIntervals = mcCycles / numprocs;
    int myloopBegin = my_rank * noIntervals + 1;
    int myloopEnd = (my_rank + 1) * noIntervals;
    if ((my_rank == numprocs - 1) && (myloopEnd < mcCycles)) myloopEnd = mcCycles;

    // Broadcast to all nodes common variables
    MPI_Bcast (&nSpins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&tempInit, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&tempFinal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&tempStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Declare new file name and add lattice size to file name
    string fileout = "Lattice";

    // Add MC cycles to filename
    //fileout.append(to_string(nSpins) + "x" + to_string(nSpins) + "_1e");
    //stringstream mcCyclesString;
    //mcCyclesString << log10(mcCycles);
    //fileout.append(mcCyclesString.str() + ".txt");

    // Without MC cycles to filename
    fileout.append(to_string(nSpins) + "x" + to_string(nSpins));
    fileout.append("_MPI.txt");

    ofile.open(fileout);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << "MC cycles:" << setw(15) << "Temperature:" << setw(15) << "Energy:" <<
             setw(15) << "Cv:" << setw(15) << "Magnetization:" << setw(15) << "Chi (X):" <<
             setw(15) << "|Magnetization|:" << setw(15) << "# Flips:" << endl;

    double timeStart = MPI_Wtime();
    double CvError, XError;
    for(int i = 0; i < numberOfLoops; i++)
    {
        // Start Metropolis algorithm with Monte Carlo sampling
        for(double temp = tempInit; temp <= tempFinal; temp += tempStep)
        {
            double acceptedFlips = 0.0;
            double newAcceptedFlips = 0.0;
            if(my_rank == 0)
            {
                cout << "Temperature: " << temp << endl;
            }
            vec expectationValues = zeros<mat>(6);
            vec totalExpectationValues = zeros<mat>(6);

            double energy = 0.0;
            double magneticMoment = 0.0;
            mat spinMatrix = zeros<mat>(nSpins,nSpins);
            initializeLattice(nSpins, spinMatrix, energy, magneticMoment);

            vec energyDifference = zeros<mat>(17);
            for(int dE = -8; dE <= 8; dE += 4)
            {
                energyDifference(dE + 8) = exp(-dE/temp);
            }

            for(double cycles = myloopBegin; cycles <= myloopEnd; cycles++)
            {
                metropolis(nSpins, spinMatrix, energy, magneticMoment, acceptedFlips, energyDifference);

                expectationValues(0) += energy;
                expectationValues(1) += energy * energy;
                expectationValues(2) += fabs(energy);
                expectationValues(3) += magneticMoment;
                expectationValues(4) += magneticMoment * magneticMoment;
                expectationValues(5) += fabs(magneticMoment);
            }

            // Find total average
            for(int i = 0; i < 6; i++)
            {
                MPI_Reduce(&expectationValues[i], &totalExpectationValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }
            MPI_Reduce(&acceptedFlips, &newAcceptedFlips, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            newAcceptedFlips /= mcCycles;

            if(my_rank == 0)
            {
                writeToFile(nSpins, mcCycles, temp, totalExpectationValues, CvError, XError, newAcceptedFlips);
                cout << "Accepted flips: " << newAcceptedFlips << endl;
            }
        }
    }

    CvError /= numberOfLoops;
    XError /= numberOfLoops;

    //cout << "Cv error: " << CvError << endl;
    //cout << "X error: " << XError << endl;

    ofile.close();

    double timeEnd = MPI_Wtime();
    double totalTime = timeEnd - timeStart;
    if (my_rank == 0)
    {
        cout << "Time = " <<  totalTime  << " seconds on number of processors: "  << numprocs  << endl;
    }

    // End MPI
    MPI_Finalize ();

    return 0;
}

// Initalize lattice with spins, magnetic moment and energy
void initializeLattice(int nSpins, mat &spinMatrix, double &energy, double &magneticMoment)
{
    // Initialize RNG, can be called by rand(gen) to get a random number between 0 and 1
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> rand(0.0, 1.0);

    // Spin state sets the ground state of the spins, 1 is all up, -1 is all down and 0 is random
    int spinState = 0;

    // Set ground state of the lattice
    for(int y = 0; y < nSpins; y++)
    {
        for(int x = 0; x < nSpins; x++)
        {
            if(spinState == 1) spinMatrix(x,y) = 1.0; // All spins are up
            if(spinState == -1) spinMatrix(x,y) = -1.0; // All spins are down
            if(spinState == 0) // Random spin directions
            {
                if(rand(gen) < 0.5) spinMatrix(x,y) = -1.0;
                if(rand(gen) >= 0.5) spinMatrix(x,y) = 1.0;
            }
        }
    }
    //cout << spinMatrix << endl;

    // Set up magnetic moment and energy
    for(int y = 0; y < nSpins; y++)
    {
        for(int x = 0; x < nSpins; x++)
        {
            magneticMoment += (double) spinMatrix(x,y);
            energy -= spinMatrix(x,y) * (spinMatrix(pbc(x,nSpins,-1),y) +
                                         spinMatrix(x,pbc(y,nSpins,-1)));
        }
    }
}

// Perform the Metropolis algorithm
void metropolis(int nSpins, mat &spinMatrix, double &energy, double &magneticMoment, double &acceptedFlips, vec &energyDifference)
{
    // Initialize RNG, can be called by rand(gen) to get a random number between 0 and 1
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> rand(0.0, 1.0);

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
                acceptedFlips++;
            }
        }
    }

    //cout << "Energy: " << energy << endl;
    //cout << "Magnetic moment: " << magneticMoment << endl;
    //cout << "MC cycles: " << mcCycles << endl;
}

void writeToFile(int nSpins, double mcCycles, double temp, vec expectationValues, double &CvError, double &XError, double acceptedFlips)
{
    // Normalization of the values
    double norm = 1.0 / mcCycles;
    double normSpins = 1.0 / nSpins / nSpins;

    // Numerical values
    double expectVal_E = expectationValues(0) * norm;// / nSpins / nSpins;
    double expectVal_E2 = expectationValues(1) * norm;
    double expectVal_Eabs = expectationValues(2) * norm;
    double expectVal_M = expectationValues(3) * norm;// / nSpins / nSpins;
    double expectVal_M2 = expectationValues(4) * norm;// / nSpins / nSpins;
    double expectVal_Mabs = expectationValues(5) * norm;

    double expectVal_Cv = (expectVal_E2 - expectVal_Eabs * expectVal_Eabs) * normSpins / (temp * temp);
    double expectVal_X = (expectVal_M2 - expectVal_Mabs * expectVal_Mabs) * normSpins / temp;

    //cout << "\nExpectation values, numerical: " << endl;
    //cout << "Energy: " << expectVal_E << endl;
    //cout << "Energy^2: " << expectationValues_E2 << endl;
    //cout << "|Energy|: " << expectVal_Eabs << endl;
    //cout << "Magnetic moment: " << expectVal_M << endl;
    //cout << "Magnetic moment^2: " << expectationValues_M2 << endl;
    //cout << "|Magnetic moment|: " << expectVal_Mabs << endl;
    //cout << "Heat capacity: " << expectVal_Cv << endl;
    //cout << "Susceptibility: " << expectVal_X << endl;

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

    //cout << "\nExpectation values, analytical:" << endl;
    //cout << "Energy: " << expectValAnalytical_E << endl;
    //cout << "Heat capacity: " << expectValAnalytical_Cv << endl;
    //cout << "Susceptibility: " << expectValAnalytical_X << endl;

    // Error between numerical and analytical
    //vec CvError(10);
    //vec XError(10);
    CvError += fabs(expectValAnalytical_Cv - expectVal_Cv) / expectValAnalytical_Cv;
    XError += fabs(expectValAnalytical_X - expectVal_X) / expectValAnalytical_X;
    //cout << "Cv error: " << CvError << endl;
    //cout << "X error: " << XError << endl;

    //ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << mcCycles;
    ofile << setw(15) << setprecision(8) << temp;
    ofile << setw(15) << setprecision(8) << expectVal_E * normSpins;
    ofile << setw(15) << setprecision(8) << expectVal_Cv;
    ofile << setw(15) << setprecision(8) << expectVal_M * normSpins;
    ofile << setw(15) << setprecision(8) << expectVal_X;
    ofile << setw(15) << setprecision(8) << expectVal_Mabs * normSpins;
    ofile << setw(15) << setprecision(8) << acceptedFlips << endl;
}
