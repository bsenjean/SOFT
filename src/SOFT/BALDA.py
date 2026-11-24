#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np

def correlation_BALDA(U,t,n,beta,dbeta_dU):
   """
   Function that generates the correlation BALDA potential and energy for one site with occupation n.
   It also generates the derivative of the correlation energy with respect to U and t, which can be useful
   to compute some particular quantities (like the double occupation).
   """

   ts = - 4.0*t*np.sin(np.pi*n/2.0)/np.pi
   ec = - 2.*t*beta*np.sin(np.pi*n/beta)/np.pi + 4.0*t*np.sin(np.pi*n/2.0)/np.pi - U*0.25*n*n
   dec_dU = dbeta_dU*(-2.*t*np.sin(np.pi*n/beta)/np.pi + 2.*t*n*np.cos(np.pi*n/beta)/beta) - n*n*0.25
   dec_dn = - 2.*t*np.cos(np.pi*n/beta) + 2.*t*np.cos(np.pi*n/2.0) - U*0.5*n
   dts_dn = - 2.*t*np.cos(np.pi*n/2.0)

   # Hole-particle symmetry:
   if n > 1.0:
       ts = - 4.0*t*np.sin(np.pi*(2.-n)/2.0)/np.pi
       ec = - 2.*t*beta*np.sin(np.pi*(2.-n)/beta)/np.pi + U*(n-1.) + 4.0*t*np.sin(np.pi*(2.-n)/2.0)/np.pi - U*n*n*0.25
       dec_dU = dbeta_dU*(-2.*t*np.sin(np.pi*(2.-n)/beta)/np.pi + 2.*t*(2.-n)*np.cos(np.pi*(2.-n)/beta)/beta) + n - 1. - n*n*0.25
       dec_dn = 2.*t*np.cos(np.pi*(2.0-n)/beta) - 2.*t*np.cos(np.pi*(2.0-n)/2.0) + U - U*0.5*n
       dts_dn = 2.*t*np.cos(np.pi*(2.0-n)/2.0)

   dec_dt = ec/t - U*dec_dU/t

   return ts, ec, dec_dU, dec_dt, dec_dn, dts_dn

# ---------- Homogeneous 2D reference: t_s(n) as a numerical table ----------

def build_ts_table_2D(t=1.0, nk=2000, n_grid=None):
    """Per-site non-interacting kinetic energy t_s(n) of the 2D square lattice,
    t_s(n) = -(4t/pi)*... (no closed form in 2D) -> computed numerically."""
    if n_grid is None:
        n_grid = np.linspace(1e-6, 2 - 1e-6, 400)
    k = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    KX, KY = np.meshgrid(k, k, indexing='ij')
    eps = -2 * t * (np.cos(KX) + np.cos(KY))          # dispersion
    eps_flat = eps.ravel()
    dA = (2 * np.pi / nk) ** 2                        # k-space volume element
    mu_grid = np.linspace(eps_flat.min(), eps_flat.max(), 4000)
    n_mu  = np.array([2 * dA * np.sum(eps_flat <= mu) for mu in mu_grid])   # spin factor 2
    ts_mu = np.array([dA * np.sum(eps_flat[eps_flat <= mu]) for mu in mu_grid])
    # invert n(mu) -> t_s(n) by interpolation
    ts_n = np.interp(n_grid, n_mu, ts_mu)
    mu_n  = np.interp(n_grid, n_mu, mu_grid)   # <-- d t_s / d n
    return n_grid, ts_n, mu_n

# ---------- Local (LDA) kinetic energy per site for a finite lattice ----------

def ts_per_site(n, n_grid, ts_n):
    """n: array of local densities n_i in [0,2] (one per site).
    Returns t_s(n_i) per site and the total T_s^LDA."""
    n = np.clip(n, n_grid[0], n_grid[-1])
    ts_sites = np.interp(n, n_grid, ts_n)
    return ts_sites, np.sum(ts_sites)

def dts_dn_per_site(n, n_grid, mu_n):
    """d t_s / d n at each site = mu_s(n_i)."""
    n = np.clip(n, n_grid[0], n_grid[-1])
    return np.interp(n, n_grid, mu_n)

def generate_potential_1D(L,U,t,occ,beta,dbeta_dU):
   """
   Function that generates the Hxc BALDA potentials and energies for the full system.
   """

   deHxc_dn = np.zeros((L),dtype=float)
   dts_dn   = np.zeros((L),dtype=float)
   eHxc     = np.zeros((L),dtype=float)
   ts       = np.zeros((L),dtype=float)

   for i in range(L):
     ts[i], ec, dec_dU, dec_dt, dec_dn, dts_dn[i] = correlation_BALDA(U,t,occ[i],beta,dbeta_dU)
     eHxc[i]     = ec + U*0.25*occ[i]*occ[i]
     deHxc_dn[i] = dec_dn + U*occ[i]*0.5

   return ts, eHxc, deHxc_dn, dts_dn

def generate_potential_2D(Lx,Ly,U,t,opdm,t_matrix,n_grid,ts_table,dts_dn_table,beta,dbeta_dU):
   """
   Function that generates the Hxc dimensional-interpolation-ansatz potentials and energies for the full system.
   """

   """
   Based on the work https://doi.org/10.1088/1361-648X/ab32ac
   Eq. 18
   The correlation functional from the 2D model can be approximated from the 1D model with a Viriel-like formula
   e(d,U,n) ~ e(d,0,n)*e(d=1,U,n)/e(d=1,0,n), where e(d=1,0,n) is the non-interacting part of the 1D model, and e(d=1,U,n) the interacting one, both can be approximated BALDA.
   e(d,0,n) is the only problem, but it should be equal to the sum over occupied orbital energies. The derivative wrt n should be the chemical potential, i.e. the HOMO orbital energy.
   With e(d,U,n) = e(d,0,n) + (U/4)n**2 + ec(d,U,n) ; we finally have:
   ec(d,U,n) = e(d,0,n)*(e(1,U,n)/e(1,0,n) - 1) - (U/4)n**2

   Actually we can simplify further:
   eHxc(d) = ts(d)*eHxc(d=1)/ts(d=1)
   
   and

   deHxc_dn(d) = dts_dn(d)*(eHxc(d=1)/ts(d=1)) + ts(d)*(deHxc_dn(d=1)*ts(d=1) - eHxc(d=1)*dts_dn(d=1))/(ts(d=1)**2) 

   where ts(d) = e(d,U=0)
   """

   L = Lx*Ly
   occ = np.diag(opdm)
   eHxc_2D     = np.zeros((L),dtype=float)
   deHxc_dn_2D = np.zeros((L),dtype=float)

   Ts_2D = np.trace(t_matrix @ opdm) # just to compare with the LDA formula?

   # LDA based on tabulated values of ts and dts_dn:
   ts_2D, Ts_LDA = ts_per_site(occ, n_grid, ts_table)
   dts_dn_2D = dts_dn_per_site(occ, n_grid, dts_dn_table)

   # comparison. Should be the same for a uniform system:
   print(f"T_s^LDA total  = {Ts_LDA:.6f}")
   print(f"T_s^2D total  = {Ts_2D:.6f}")

   # BALDA, for the 1D model:
   ts_1D, eHxc_1D, deHxc_dn_1D, dts_dn_1D = generate_potential_1D(L,U,t,occ,beta,dbeta_dU)

   # 2D model: dimensional-interpolation ansatz:
   for i in range(L):
     eHxc_2D[i] = ts_2D[i]*eHxc_1D[i]/ts_1D[i]
     deHxc_dn_2D[i] = dts_dn_2D[i]*(eHxc_1D[i]/ts_1D[i]) + ts_2D[i]*(deHxc_dn_1D[i]*ts_1D[i] - eHxc_1D[i]*dts_dn_1D[i])/(ts_1D[i]**2)

   return ts_2D, eHxc_2D, deHxc_dn_2D
