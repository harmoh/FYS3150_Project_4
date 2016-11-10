# Runs main.cpp and plots the generated .txt files

import os
from math import *
import numpy as np
from matplotlib import pyplot as plt

os.system('c++ main.cpp -o main.o -O3 -I /usr/local/Cellar/armadillo/7.400.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack')

def read(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    x0 = []; x1 = []; x2 = []; x3 = []; x4 = []; x5 = []; x6 = []; x7 = [];
    # Read lines except for the first one:
    lines = infile.readlines()[1:]
    for line in lines:
        words = line.split()
        x0.append(float(words[0]))
        x1.append(float(words[1])) # Temperature
        x2.append(float(words[2]))
        x3.append(float(words[3]))
        x4.append(float(words[4]))
        x5.append(float(words[5]))
        x6.append(float(words[6]))
        x7.append(float(words[7])) # Number of flips
    infile.close()
    return x1, x2, x3, x4, x5, x6, x7

run = './main.o 20 10000 1.0 2.5 0.1 1';
os.system(run) # Argument for number of spins, MC cycles, initial and final temperature, tempurate step and number of loops.

# Fetching data by a call on read_x_u_v for three different n:
x1, x2, x3, x4, x5, x6, x7 = read('Lattice20x20.txt')

plt.xlabel('Temperature')
plt.ylabel('Number of accepted flips')
#plt.xscale('log', nonposy='clip')
plt.yscale('log', nonposy='clip')
plt.rcParams.update({'font.size': 10})
#plt.axis([10, 1000000, -1.99, -2])
axes = plt.gca()
axes.set_xlim([1.0, 2.4])
#axes.set_ylim([-2,-1.99])
plt.plot(x1, x7, linewidth = 1.0, label = '# Accepted flips')
#plt.plot(x1, x3, linewidth = 1.0, label = 'Magnetic moment')
#plt.plot(x1, x4, label = '$\chi$, numerical')
#plt.plot(x1, analytical_Cv, label = '$C_V$, analytical')
#plt.plot(x1, analytical_X, label = '$\chi$, analytical')
plt.legend(loc='upper right',fancybox='True')
plt.grid()
plt.savefig('Lattice20x20_temp_.eps', format = 'eps', dpi = 1000, bbox_inches='tight')
#plt.show();
