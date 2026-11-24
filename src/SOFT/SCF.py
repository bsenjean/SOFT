#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import numpy as np
from scipy.linalg import eigh
from scipy import optimize
import random

sys.path.insert(0, os.getcwd())
from .BALDA import build_ts_table_2D
from .BALDA import generate_potential_1D
from .BALDA import generate_potential_2D

def generate_hamiltonian(Lx,Ly,N,t,v,periodic=True):
   """
   Function that generates the non-interacting 1D or 2D Hubbard Hamiltonian with potential v = [v_00, v_01, v_02, ... , v_10, v_11, ..., v_(Lx-1,Ly-1)]
   Each site will be labelled by a single index p=x+Lx(y−1),x=1..Lx,y=1..Ly (or p=x+Lx*y if we start at 0)
   Periodic boundary condition: add -t to ±(L−Lx+1)-type positions, i.e. between (Lx,y) and (1,y), and between (x,Ly) and (x,1)
   """

   L = Lx*Ly
   H = np.zeros((L,L),dtype=float) # could be H = np.diag(v.copy())

   for y in range(Ly):
     for x in range(Lx):
       p = x+Lx*y
       H[p,p] = v[p]
       if x + 1 < Lx: # neighbor along +x
         H[p,p+1] = H[p+1,p] = -t
       elif periodic and Lx > 1:
         H[p,y*Lx] = H[y*Lx,p] = -t
       if y + 1 < Ly: # neighbor along +y
         H[p,p+Lx] = H[p+Lx,p] = -t
       elif periodic and Ly > 1:
         H[p,x] = H[x,p] = -t

   # For 1D, sometimes we need antiperiodic condition, I add it here:
   if Ly==1 and ((N//2) % 2 == 0): #antiperiodic condition
     H[L-1,0] = H[0,L-1] = +t

   return H

def solve_hamiltonian(L,N,H):
   """
   Function to diagonalize the Hamiltonian and extract the one-particle reduced density matrix and the orbital energies.
   """

   orb_energy,C = eigh(H)

   # One-particle reduced density matrix:
   onepdm = np.zeros((L,L),dtype=float)
   occ = np.zeros((L),dtype=float)
   for i in range(L):
    for j in range(L):
     for k in range(N//2):
      onepdm[i,j]+=2*C[i,k]*C[j,k] # factor 2 because we work with spatial orbitals
    occ[i] = onepdm[i,i] 

   return onepdm,occ,orb_energy,C

def self_consistency(Lx,Ly,Nelec,U,t,pot,beta,dbeta_dU,MAXIT,path_results_folder,output_file,mix_cst=0.4):
   """
   Function to get the self-consistent energy and occupations.
   """

   Nsites = Lx*Ly

   # Start with the opdm from the non-interacting Hamiltonian
   t_matrix = generate_hamiltonian(Lx,Ly,Nelec,t,np.zeros(Nsites),periodic=True)
   H0 = generate_hamiltonian(Lx,Ly,Nelec,t,pot)
   opdm, occ, orb_energy, C = solve_hamiltonian(Nsites,Nelec,H0)

   print("Fermi gap of H0:")
   print(orb_energy[Nelec//2] - orb_energy[Nelec//2 - 1])

   with open(output_file,'w') as f:
     f.write("*"*10 + " Initialization " + "*"*10 + "\n")
     f.write("""
 Lx= {}
 Ly= {}
 N = {}
 t = {}
 U = {}
 v = {}\n""".format(Lx,Ly,Nelec,t,U,pot))

     f.write(" Initial trial sites occupation:\n{}\n\n".format(occ))

   Etot    = 0
   it      = 0
   delta_E = 1
   convergence = True

   # If 2D, tabulate the non-interacting kinetic density-functional
   if Ly > 1: n_grid, ts_table, dts_dn_table = build_ts_table_2D(t=t)

   while (delta_E > 1e-4 or normocc > 1e-4) and convergence: #or mix_cst < 1) and convergence:


     # Compute the vHxc contributions 

     if Ly == 1:
       deHxc_dn = generate_potential_1D(Lx,U,t,occ,beta,dbeta_dU)[2]
     else:
       deHxc_dn = generate_potential_2D(Lx,Ly,U,t,opdm,t_matrix,n_grid,ts_table,dts_dn_table,beta,dbeta_dU)[2]
     KSpot = np.add(pot,deHxc_dn)

     # Construct and solve the KS Hamiltonian
     h_SOFT = generate_hamiltonian(Lx,Ly,Nelec,t,KSpot)
     opdm, occ_new, KS_energies, KS_orbs = solve_hamiltonian(Nsites,Nelec,h_SOFT)

     # compute the eHxc contribution with the new density

     # BALDA in 1D:
     if Ly == 1: 
       eHxc = generate_potential_1D(Nsites,U,t,occ_new,beta,dbeta_dU)[1]
     else:
       eHxc = generate_potential_2D(Lx,Ly,U,t,opdm,t_matrix,n_grid,ts_table,dts_dn_table,beta,dbeta_dU)[1]

     # Compute the total energy
     sum_occ_eKS  = 0.
     sum_eHxc     = 0.
     sum_deHxc_dn = 0.
     for i in range(Nelec//2): sum_occ_eKS += 2*KS_energies[i]
     for i in range(Nsites): sum_eHxc += eHxc[i]
     for i in range(Nsites): sum_deHxc_dn += deHxc_dn[i]*occ_new[i]
     Etot_new = sum_occ_eKS + sum_eHxc - sum_deHxc_dn 

     delta_E = abs(Etot_new - Etot)
     Etot    = Etot_new
     normocc = np.linalg.norm(occ - occ_new)

     #if (delta_E <= 1e-4 and normocc <= 1e-4): mix_cst += 0.2
     occ = (1 - mix_cst)*occ + mix_cst*occ_new

     with open(output_file,'a') as f:
       f.write("*"*10 + " ITERATION {:3d} ".format(it) + "*"*10 + "\n")
       f.write("mix constant     : {:4.2f}".format(mix_cst) + "\n")
       f.write("Energy (hartree) : {:16.8f}".format(Etot) + "\n")
       f.write("Occupied KS      : {}".format(KS_energies[:Nelec//2]) + "\n")
       f.write("New    occ       : {}".format(occ_new) + "\n")
       f.write("Damped occ       : {}".format(occ) + "\n")
       f.write("Delta Energy     : {:16.8f}".format(delta_E) + "\n")
       f.write("Norm Delta_occ   : {:16.8f}\n".format(normocc) + "\n")

     it+=1

     if it > MAXIT: convergence=False

   return convergence, Etot, occ, it, KS_orbs
    
def run_SOFT_Hubbard(Lx,Ly,Nelec,U,t,pot,MAXIT,code_directory,output_file,mix_cst=0.4):

   # Compute beta and dbeta_dU for the approximate BALDA potential
   subprocess.check_call("echo " + str(U) + " " + str(t) + " | beta_and_derivatives",shell=True, cwd = code_directory)
   with open(code_directory + "beta_dbetadU.dat","r") as f:
      line = f.read()
      beta = float(line.split()[0])
      dbeta_dU = float(line.split()[1])
      f.close()

   # Perform the self-consistent algorithm
   conv, Etot, occ, it, KS_orbs = self_consistency(Lx,Ly,Nelec,U,t,pot,beta,dbeta_dU,MAXIT,code_directory+"/results/",output_file,mix_cst=mix_cst)

   # Check convergence.
   if (conv):
     with open(output_file,'a') as f:
        f.write("*"*10 + " SUCCESS " + "*"*10 + "\n")
        f.write("Iteration        : {:16d}".format(it) + "\n")
        f.write("Energy (hartree) : {:16.8f}".format(Etot) + "\n")
        f.write("Occupations      : {}".format(occ) + "\n")
   else:
     with open(output_file,'a') as f:
        f.write("*"*10 + " FAILURE " + "*"*10 + "\n")
        f.write("Iteration > {}, NO_CONVERGENCE_REACHED".format(MAXIT) + "\n")

   return KS_orbs, Etot, occ
