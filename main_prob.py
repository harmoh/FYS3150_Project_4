# Runs main.cpp and plots the generated .txt files

import os
from math import *
import numpy as np
from matplotlib import pyplot as plt

os.system('mpic++ main.cpp -o main.o -O3 -I /usr/local/Cellar/armadillo/7.400.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack')

def read(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    x0 = []; x1 = [];
    # Read lines except for the first one:
    lines = infile.readlines()[1:]
    for line in lines:
        words = line.split()
        x0.append(float(words[0])) # Energy
    infile.close()
    return x0

run = 'mpirun -n 1 ./main.o 20 10000000 1.0 1.0 0.1 1';
os.system(run) # Argument for number of spins, MC cycles, initial and final temperature, tempurate step and number of loops.

# Fetching data by a call on read_x_u_v for three different n:
x0 = read('Probability.txt')

plt.xlabel('Energy')
plt.ylabel('Number of values (normalized)')
plt.rcParams.update({'font.size': 10})
plt.hist([x0], normed = 1, label = '$T = 1.0$') # Change to correct temperature
plt.legend(loc='upper right',fancybox='True')
plt.grid()
plt.savefig('Lattice20x20_prob_t1.0_1e7_.eps', format = 'eps', dpi = 1000, bbox_inches='tight') # Change to appropriate name
#plt.show();
