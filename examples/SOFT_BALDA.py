#!/usr/bin/python 
import numpy as np
import sys
import os
import SOFT
import random

working_directory = os.getenv('SOFT_DIR') + "/examples/"

n_sites     = 100
n_elec      = 30
U           = 4.0
t           = 1.0
v_choice    = ["uniform","random","ABAB","power","decreasing"][3]
SCF_maxiter = 100

v = np.full((n_sites),1.*n_elec/(1.*n_sites)) # uniform v_choice
if v_choice == "uniform":
  pass
elif v_choice == "random":
  for i in range(n_sites): v[i] = random.uniform(-1,1)
elif v_choice == "ABAB":
  for i in range(n_sites//2):
    v[2*i]   = -1.
    v[2*i+1] = +1.
elif v_choice == "power": # See GAO XIANLONG et al. PHYSICAL REVIEW B 73, 165120 (2006)
  l = 2
  const_V = 0.006
  for i in range(n_sites):
    v[i] = const_V*(i - n_sites/2)**l
elif v_choice == "decreasing":
  for i in reversed(range(n_sites)):
     v[i] = 0.1*i
else:
  sys.exit("The v_choice is not defined. Program terminated.")

output_file = working_directory + "/results/L{}_N{}_U{}_t{}_{}.dat".format(n_sites,n_elec,U,t,v_choice)

KS_orbs, E_DFT, density_exact = SOFT.run_SOFT_Hubbard(n_sites,n_elec,U,t,v,SCF_maxiter,working_directory,output_file)
