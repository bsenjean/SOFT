#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np

def correlation_BALDA(U,t,n,beta,dbeta_dU):
   """
   Function that generates the correlation BALDA potential and energy for one site with occupation n.
   It also generates the derivative of the correlation energy with respect to U and t, which can be useful
   to compute some particular quantities (like the double occupation).
   """

   ec = - 2.*t*beta*np.sin(np.pi*n/beta)/np.pi + 4.0*t*np.sin(np.pi*n/2.0)/np.pi - U*0.25*n*n
   dec_dU = dbeta_dU*(-2.*t*np.sin(np.pi*n/beta)/np.pi + 2.*t*n*np.cos(np.pi*n/beta)/beta) - n*n*0.25
   dec_dn = - 2.*t*np.cos(np.pi*n/beta) + 2.*t*np.cos(np.pi*n/2.0) - U*0.5*n

   # Hole-particle symmetry:
   if n > 1.0:
       ec = - 2.*t*beta*np.sin(np.pi*(2.-n)/beta)/np.pi + U*(n-1.) + 4.0*t*np.sin(np.pi*(2.-n)/2.0)/np.pi - U*n*n*0.25
       dec_dU = dbeta_dU*(-2.*t*np.sin(np.pi*(2.-n)/beta)/np.pi + 2.*t*(2.-n)*np.cos(np.pi*(2.-n)/beta)/beta) + n - 1. - n*n*0.25
       dec_dn = 2.*t*np.cos(np.pi*(2.0-n)/beta) - 2.*t*np.cos(np.pi*(2.0-n)/2.0) + U - U*0.5*n

   dec_dt = ec/t - U*dec_dU/t

   return ec, dec_dU, dec_dt, dec_dn


def generate_potential(L,U,t,occ,beta,dbeta_dU):
   """
   Function that generates the Hxc BALDA potentials and energies for the full system.
   """

   deHxc_dn = np.zeros((L),dtype=float)
   deHxc_dt = np.zeros((L),dtype=float)
   deHxc_dU = np.zeros((L),dtype=float)
   eHxc     = np.zeros((L),dtype=float)

   for i in range(L):
     ec, dec_dU, dec_dt, dec_dn = correlation_BALDA(U,t,occ[i],beta,dbeta_dU)
     eHxc[i]     = ec + U*0.25*occ[i]*occ[i]
     deHxc_dn[i] = dec_dn + U*occ[i]*0.5

   return eHxc, deHxc_dn
