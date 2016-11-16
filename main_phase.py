# Runs main.cpp and plots the generated .txt files

import os
from math import *
import numpy as np
from matplotlib import pyplot as plt

os.system('mpic++ main.cpp -o main.o -O3 -I /usr/local/Cellar/armadillo/7.400.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack')

def read(filename):
    infile = open(filename, 'r')
    # Elements to be read in file:
    x1 = []; x2 = []; x3 = []; x4 = []; x5 = [];
        
    # Read lines except for the first one:
    lines = infile.readlines()[1:]
    for line in lines:
        words = line.split()
        x1.append(float(words[1])) # Temp
        x2.append(float(words[2])) # <E>
        x3.append(float(words[3])) # <|M|>
        x4.append(float(words[4])) # Cv
        x5.append(float(words[5])) # X
    infile.close()
    return x1, x2, x3, x4, x5

x40_1 = []; x40_2 = []; x40_3 = []; x40_4 = []; x40_5 = []; x60_1 = []; x60_2 = []; x60_3 = []; x60_4 = []; x60_5 = []; x100_1 = []; x100_2 = []; x100_3 = []; x100_4 = []; x100_5 = []; x140_1 = []; x140_2 = []; x140_3 = []; x140_4 = []; x140_5 = [];

mcCyclesExponential = 6 # Written as the exponential
mcCycles = 10**mcCyclesExponential;
temp_init = 2.18;
temp_final = 2.32;
temp_step = 0.01;

# Argument for number of spins, MC cycles, initial and final temperature, tempurate step and number of loops.
# L = 40
run = 'mpirun -n 4 ./main.o 40 ' + str(mcCycles) + ' ' + str(temp_init) + ' ' + str(temp_final) + ' ' + str(temp_step) + ' 1';
os.system(run)
# Fetching data by a call on read()
x40_1, x40_2, x40_3, x40_4, x40_5 = read('Phase40.txt')

# L = 60
run = 'mpirun -n 4 ./main.o 60 ' + str(mcCycles) + ' ' + str(temp_init) + ' ' + str(temp_final) + ' ' + str(temp_step) + ' 1';
os.system(run)
x60_1, x60_2, x60_3, x60_4, x60_5 = read('Phase60.txt')

# L = 100
run = 'mpirun -n 4 ./main.o 100 ' + str(mcCycles) + ' ' + str(temp_init) + ' ' + str(temp_final) + ' ' + str(temp_step) + ' 1';
os.system(run)
# Fetching data by a call on read()
x100_1, x100_2, x100_3, x100_4, x100_5 = read('Phase100.txt')

# L = 140
run = 'mpirun -n 4 ./main.o 140 ' + str(mcCycles) + ' ' + str(temp_init) + ' ' + str(temp_final) + ' ' + str(temp_step) + ' 1';
os.system(run)
x140_1, x140_2, x140_3, x140_4, x140_5 = read('Phase140.txt')

plt.figure()
plt.xlabel('Temperature')
plt.ylabel(r'$\langle  E \rangle / L^2$')
plt.rcParams.update({'font.size': 10})
axes = plt.gca()
axes.set_xlim([temp_init, temp_final+0.00001])
#axes.set_ylim([-2,-1.99])
plt.plot(x40_1, x40_2, linewidth = 1.0, label = '$L = 40$')
plt.plot(x60_1, x60_2, linewidth = 1.0, label = '$L = 60$')
plt.plot(x100_1, x100_2, linewidth = 1.0, label = '$L = 100$')
plt.plot(x140_1, x140_2, linewidth = 1.0, label = '$L = 140$')
plt.legend(loc='upper left',fancybox='True')
plt.grid()
plt.savefig('Phase_1e' + str(mcCyclesExponential) + '_energy_.eps', format = 'eps', dpi = 1000, bbox_inches='tight')

plt.figure()
plt.xlabel('Temperature')
plt.ylabel(r'$\langle  M \rangle / L^2$')
plt.rcParams.update({'font.size': 10})
axes = plt.gca()
axes.set_xlim([temp_init, temp_final+0.00001])
#axes.set_ylim([-2,-1.99])
plt.plot(x40_1, x40_3, linewidth = 1.0, label = '$L = 40$')
plt.plot(x60_1, x60_3, linewidth = 1.0, label = '$L = 60$')
plt.plot(x100_1, x100_3, linewidth = 1.0, label = '$L = 100$')
plt.plot(x140_1, x140_3, linewidth = 1.0, label = '$L = 140$')
plt.legend(loc='upper right',fancybox='True')
plt.grid()
plt.savefig('Phase_1e' + str(mcCyclesExponential) + '_mag_.eps', format = 'eps', dpi = 1000, bbox_inches='tight')

plt.figure()
plt.xlabel('Temperature')
plt.ylabel(r'$C_V / L^2$')
plt.rcParams.update({'font.size': 10})
axes = plt.gca()
axes.set_xlim([temp_init, temp_final+0.00001])
#axes.set_ylim([-2,-1.99])
plt.plot(x40_1, x40_4, linewidth = 1.0, label = '$L = 40$')
plt.plot(x60_1, x60_4, linewidth = 1.0, label = '$L = 60$')
plt.plot(x100_1, x100_4, linewidth = 1.0, label = '$L = 100$')
plt.plot(x140_1, x140_4, linewidth = 1.0, label = '$L = 140$')
plt.legend(loc='upper left',fancybox='True')
plt.grid()
plt.savefig('Phase_1e' + str(mcCyclesExponential) + '_Cv_.eps', format = 'eps', dpi = 1000, bbox_inches='tight')

plt.figure()
plt.xlabel('Temperature')
plt.ylabel(r'$\chi / L^2$')
plt.rcParams.update({'font.size': 10})
axes = plt.gca()
axes.set_xlim([temp_init, temp_final+0.00001])
#axes.set_ylim([-2,-1.99])
plt.plot(x40_1, x40_5, linewidth = 1.0, label = '$L = 40$')
plt.plot(x60_1, x60_5, linewidth = 1.0, label = '$L = 60$')
plt.plot(x100_1, x100_5, linewidth = 1.0, label = '$L = 100$')
plt.plot(x140_1, x140_5, linewidth = 1.0, label = '$L = 140$')
plt.legend(loc='upper left',fancybox='True')
plt.grid()
plt.savefig('Phase_1e' + str(mcCyclesExponential) + '_X_.eps', format = 'eps', dpi = 1000, bbox_inches='tight')



