# Runs main.cpp and plots the generated .txt files

import os
from math import *
import numpy as np
from matplotlib import pyplot as plt

os.system('mpic++ main.cpp -o main.o -O3 -I /usr/local/Cellar/armadillo/7.400.2/include -DARMA_DONT_USE_WRAPPER -lblas -llapack')

mcCyclesExponential = 3 # Written as the exponential
mcCycles = 10**mcCyclesExponential;
temp_init = 2.18;
temp_final = 2.32;
temp_step = 0.01;

# L = 140
run = 'mpirun -n 1 ./main.o 140 ' + str(mcCycles) + ' ' + str(temp_init) + ' ' + str(temp_final) + ' ' + str(temp_step) + ' 1';
os.system(run)
