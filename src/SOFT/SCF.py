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
from .BALDA import generate_potential

def generate_hamiltonian(L,N,t,v):
   """
   Function that generates the non-interacting 1D Hubbard Hamiltonian with potential v.
   """

   H = np.zeros((L,L),dtype=float)

   for i in range(L-1):
    H[i,i+1] = H[i+1,i] = -t
    H[i,i] = v[i]

   if ((N//2) % 2 == 1): #periodic condition
    H[L-1,0] = H[0,L-1] = -t
   elif ((N//2) % 2 == 0): #antiperiodic condition
    H[L-1,0] = H[0,L-1] = +t
   H[L-1,L-1] = v[L-1]

   return H

def solve_KShamiltonian(L,N,H):
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

def self_consistency(Nsites,Nelec,U,t,pot,beta,dbeta_dU,MAXIT,path_results_folder,output_file,mix_cst=0.4):
   """
   Function to get the self-consistent energy and occupations.
   """

   # Start with uniform density
   occ = np.full((Nsites),1.*Nelec/(1.*Nsites))
   with open(output_file,'a') as f:
     f.write("*"*10 + " Initialization " + "*"*10 + "\n")
     f.write("""
 L = {}
 N = {}
 t = {}
 U = {}
 v = {}\n""".format(Nsites,Nelec,t,U,pot))

     f.write(" Initial trial sites occupation:\n{}\n\n".format(occ))

   Etot    = 0
   it      = 0
   delta_E = 1
   convergence = True
   while (delta_E > 1e-4 or normocc > 1e-4) and convergence: #or mix_cst < 1) and convergence:


     # Compute the vHxc contributions 
     deHxc_dn = generate_potential(Nsites,U,t,occ,beta,dbeta_dU)[1]
     KSpot = np.add(pot,deHxc_dn)

     # Construct and solve the KS Hamiltonian
     h_SOFT = generate_hamiltonian(Nsites,Nelec,t,KSpot)
     onepdm_SOFT, occ_new, KS_energies, KS_orbs = solve_KShamiltonian(Nsites,Nelec,h_SOFT)

     # compute the eHxc contribution with the new density
     eHxc = generate_potential(Nsites,U,t,occ_new,beta,dbeta_dU)[0]

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
    
def run_SOFT_Hubbard(Nsites,Nelec,U,t,pot,MAXIT,code_directory,output_file):

   # Compute beta and dbeta_dU for the approximate BALDA potential
   subprocess.check_call("echo " + str(U) + " " + str(t) + " | beta_and_derivatives",shell=True, cwd = code_directory)
   with open(code_directory + "beta_dbetadU.dat","r") as f:
      line = f.read()
      beta = float(line.split()[0])
      dbeta_dU = float(line.split()[1])
      f.close()

   # Perform the self-consistent algorithm
   conv, Etot, occ, it, KS_orbs = self_consistency(Nsites,Nelec,U,t,pot,beta,dbeta_dU,MAXIT,code_directory+"/results/",output_file)

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
