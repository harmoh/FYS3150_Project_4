# Runs main.cpp and plots the generated .txt files

import os
from math import *
import numpy as np
from matplotlib import pyplot as plt

os.system('mpic++ main.cpp -o main.o -O3 -I /usr/local/Cellar/armadillo/7.400.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack')

def read(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    x1 = []; x2 = []; x3 = []; x4 = [];
        
    # Read lines except for the first one:
    lines = infile.readlines()[1:]
    for line in lines:
        words = line.split()
        x1.append(float(words[0]))
        x2.append(float(words[2]))
        x3.append(float(words[6]))
        #x4.append(float(words[3])) # Cv
        x4.append(float(words[5])) # X
    infile.close()
    return x1, x2, x3, x4

x1 = []; x2 = []; x3 = []; x4 = []; analytical_Cv = []; analytical_X = [];
nSpins = 20;
maxMCcycles = 7; # Written as the exponential

exponentialStepSize = (i*10**exp for exp in range(2, maxMCcycles) for i in range(1, 10))
for mcCycles in exponentialStepSize:
    run = 'mpirun -n 4 ./main.o ' + str(nSpins) + ' ' + str(mcCycles) + ' 1.0 1.0 0.01 1';
    os.system(run) # Argument for number of spins, MC cycles, initial and final temperature, tempurate step and number of loops.
    # Fetching data by a call on read_x_u_v for three different n:
    x1_temp, x2_temp, x3_temp, x4_temp = read('Lattice' + str(nSpins) + 'x' + str(nSpins) + '_MPI.txt')
    #x1_temp, x2_temp, x3_temp, x4_temp = read('Lattice20x20_1e' + str(log10(mcCycles)).split('.')[0] + '.txt');
    x1.append(x1_temp);
    x2.append(x2_temp);
    x3.append(x3_temp);
    x4.append(x4_temp);
    analytical_Cv.append(0.0320823);
    analytical_X.append(0.00401074);

plt.xlabel('Monte Carlo cycles')
plt.ylabel('Normalized by the number of spins')
plt.xscale('log', nonposy='clip')
plt.rcParams.update({'font.size': 10})
#plt.axis([10, 1000000, -1.99, -2])
axes = plt.gca()
#axes.set_xlim([xmin,xmax])
#axes.set_ylim([-2,-1.99])
plt.plot(x1, x2, linewidth = 1.0, label = 'Energy')
#plt.plot(x1, x3, linewidth = 1.0, label = 'Magnetic moment')
#plt.plot(x1, x4, label = '$\chi$, numerical')
#plt.plot(x1, analytical_Cv, label = '$C_V$, analytical')
#plt.plot(x1, analytical_X, label = '$\chi$, analytical')
plt.legend(loc='upper right',fancybox='True')
plt.grid()
plt.savefig('Lattice' + str(nSpins) + 'x' + str(nSpins) + '_mag_rand_temp=1.0_.eps', format = 'eps', dpi = 1000, bbox_inches='tight')
#plt.show();
