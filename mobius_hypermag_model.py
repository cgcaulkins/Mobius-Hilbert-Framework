#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# 1) FORCE HEADLESS BACKEND *BEFORE* pyplot
import matplotlib
matplotlib.use("Agg")  # render to files; no GUI required

import numpy as np
import matplotlib.pyplot as plt
import os

# ---------- Linear algebra helpers ----------
I2 = np.eye(2, dtype=complex)
X  = np.array([[0,1],[1,0]], dtype=complex)
Y  = np.array([[0,-1j],[1j,0]], dtype=complex)
Z  = np.array([[1,0],[0,-1]], dtype=complex)

def kron3(a,b,c):
    return np.kron(np.kron(a,b),c)

def dagger(A): return A.conj().T

def partial_trace(rho, keep, dims=(2,2,2)):
    keep = tuple(sorted(keep))
    dimS, dimH, dimR = dims
    rho = rho.reshape(dimS, dimH, dimR, dimS, dimH, dimR)
    idx_map = {0:dimS, 1:dimH, 2:dimR}
    keep_dims = [idx_map[i] for i in keep]
    out = np.zeros((np.prod(keep_dims), np.prod(keep_dims)), dtype=complex)
    for iS in range(dimS):
        for iH in range(dimH):
            for iR in range(dimR):
                for jS in range(dimS):
                    for jH in range(dimH):
                        for jR in range(dimR):
                            left=[]; right=[]
                            if 0 in keep: left.append(iS); right.append(jS)
                            if 1 in keep: left.append(iH); right.append(jH)
                            if 2 in keep: left.append(iR); right.append(jR)
                            li = np.ravel_multi_index(tuple(left),  tuple(keep_dims))
                            rj = np.ravel_multi_index(tuple(right), tuple(keep_dims))
                            traced_ok = True
                            if 0 not in keep and iS!=jS: traced_ok=False
                            if 1 not in keep and iH!=jH: traced_ok=False
                            if 2 not in keep and iR!=jR: traced_ok=False
                            if traced_ok:
                                out[li,rj] += rho[iS,iH,iR,jS,jH,jR]
    return out

def von_neumann_entropy(rho, eps=1e-12):
    vals = np.linalg.eigvalsh((rho + rho.conj().T)/2)
    vals = np.clip(vals.real, 0, 1)
    vals = vals[vals>eps]
    return float(-np.sum(vals*np.log2(vals)))

def mutual_info(rho, A, B):
    rhoA  = partial_trace(rho, A)
    rhoB  = partial_trace(rho, B)
    rhoAB = partial_trace(rho, tuple(sorted(set(A)|set(B))))
    return von_neumann_entropy(rhoA) + von_neumann_entropy(rhoB) - von_neumann_entropy(rhoAB)

# ---------- Channels / Unitaries ----------
def scipy_expm(A):
    # expm via eigendecomposition (no SciPy needed)
    vals, vecs = np.linalg.eig(A)
    return (vecs @ np.diag(np.exp(vals)) @ np.linalg.inv(vecs))

def unitary_rotation(axis, angle):
    n = np.array(axis, dtype=float); n = n / (np.linalg.norm(n)+1e-15)
    H = n[0]*X + n[1]*Y + n[2]*Z
    return scipy_expm(-1j * 0.5 * angle * H)

def iSWAP_mixer(theta):
    Hmix = np.kron(X,X) + np.kron(Y,Y)
    return scipy_expm(-1j * 0.5 * theta * Hmix)

def dephase_on_R(rho, p):
    K0 = np.sqrt(1-p) * I2
    K1 = np.sqrt(p)   * np.array([[1,0],[0,0]], dtype=complex)
    K2 = np.sqrt(p)   * np.array([[0,0],[0,1]], dtype=complex)
    out = np.zeros_like(rho, dtype=complex)
    for K in (K0,K1,K2):
        U = kron3(I2, I2, K)  # act on R
        out += U @ rho @ dagger(U)
    return out

# ---------- State prep ----------
def initial_state(mu, theta):
    a = np.sqrt(max(0.0, 1.0-mu))
    b = np.sqrt(max(0.0, mu))
    v = np.zeros(8, dtype=complex)
    ph = np.exp(1j*theta)
    def idx(s,h,r): return (s<<2) + (h<<1) + r
    v[idx(0,0,0)] += a/np.sqrt(2)
    v[idx(0,1,1)] += a*ph/np.sqrt(2)
    v[idx(1,0,1)] += b/np.sqrt(2)
    v[idx(1,1,0)] += b*ph/np.sqrt(2)
    v = v / np.linalg.norm(v)
    return v.reshape(8,1)

def density_from_state(psi):
    return psi @ psi.conj().T

# ---------- Evolution ----------
def evolve(rho, muB=0.0, n_axis=(0,0,1), theta_H=0.0, p_dephase=0.0):
    UB = unitary_rotation(n_axis, angle=muB)     # magnetic rotation on R
    UR = kron3(I2, I2, UB)
    rho1 = UR @ rho @ dagger(UR)

    UM_RH = iSWAP_mixer(theta_H)                 # Möbius/hidden transfer on (R,H)
    U_total = np.kron(I2, UM_RH)                 # (I_S ⊗ U_M_RH)
    rho2 = U_total @ rho1 @ dagger(U_total)

    rho3 = dephase_on_R(rho2, p_dephase)         # dephase on R
    return rho3

# ---------- Plots ----------
def sweep_muB_plot(mu=0.5, theta=0.0, theta_H=np.pi/4, p=0.2, savepath="MI_RH_RS_vs_muB.png"):
    muBs = np.linspace(0, np.pi, 41)
    IRH, IRS = [], []
    psi = initial_state(mu, theta)
    rho0 = density_from_state(psi)
    for muB in muBs:
        rho = evolve(rho0, muB=muB, n_axis=(0,0,1), theta_H=theta_H, p_dephase=p)
        IRH.append( mutual_info(rho, (2,), (1,)) )
        IRS.append( mutual_info(rho, (2,), (0,)) )
    plt.figure(figsize=(7,4.5))
    plt.plot(muBs, IRH, label='I(R:H)')
    plt.plot(muBs, IRS, label='I(R:S)')
    plt.xlabel(r'$\mu_B$ (magnetic rotation angle)')
    plt.ylabel('Mutual Information (bits)')
    plt.title('Hypermagnetic Organization: I(R:H) vs I(R:S) across $\mu_B$')
    plt.legend(); plt.tight_layout()
    plt.savefig(savepath, dpi=160); plt.close()

def sweep_p_plot(mu=0.5, theta=0.0, theta_H=np.pi/4, muB=1.0, savepath="MI_RH_RS_vs_p_for_muB.png"):
    ps = np.linspace(0, 0.8, 41)
    IRH, IRS = [], []
    psi = initial_state(mu, theta)
    rho0 = density_from_state(psi)
    for p in ps:
        rho = evolve(rho0, muB=muB, n_axis=(0,0,1), theta_H=theta_H, p_dephase=p)
        IRH.append( mutual_info(rho, (2,), (1,)) )
        IRS.append( mutual_info(rho, (2,), (0,)) )
    plt.figure(figsize=(7,4.5))
    plt.plot(ps, IRH, label='I(R:H)')
    plt.plot(ps, IRS, label='I(R:S)')
    plt.xlabel('Dephasing on R (p)')
    plt.ylabel('Mutual Information (bits)')
    plt.title(f'Dephasing robustness at μ_B={muB:.2f}: I(R:H) vs I(R:S)')
    plt.legend(); plt.tight_layout()
    plt.savefig(savepath, dpi=160); plt.close()

def page_like_entropy(mu=0.5, theta=0.0, theta_H=np.pi/4, muB=1.0, p=0.2, steps=30, savepath="entropy_R_page_like.png"):
    psi = initial_state(mu, theta)
    rho = density_from_state(psi)
    entR = []
    for k in range(steps):
        rho = evolve(rho, muB=muB, n_axis=(0,0,1), theta_H=theta_H, p_dephase=p)
        rhoR = partial_trace(rho, (2,))
        entR.append(von_neumann_entropy(rhoR))
    plt.figure(figsize=(7,4.5))
    plt.plot(range(steps), entR, marker='o', ms=3)
    plt.xlabel('Step (toy evaporation iteration)')
    plt.ylabel('S(R) (bits)')
    plt.title('Toy Page-like evolution of S(R) under hypermagnetic + hidden transfer')
    plt.tight_layout()
    plt.savefig(savepath, dpi=160); plt.close()

# ---------- Main ----------
if __name__ == "__main__":
    MU_DEFAULT     = 0.5
    THETA_DEFAULT  = 0.0
    THETA_HIDDEN   = np.pi/4
    P_DEPH         = 0.2
    MU_B_DEFAULT   = 1.0

    # Generate plots
    sweep_muB_plot(mu=MU_DEFAULT, theta=THETA_DEFAULT,
                   theta_H=THETA_HIDDEN, p=P_DEPH,
                   savepath="MI_RH_RS_vs_muB.png")

    sweep_p_plot(mu=MU_DEFAULT, theta=THETA_DEFAULT,
                 theta_H=THETA_HIDDEN, muB=MU_B_DEFAULT,
                 savepath="MI_RH_RS_vs_p_for_muB.png")

    page_like_entropy(mu=MU_DEFAULT, theta=THETA_DEFAULT,
                      theta_H=THETA_HIDDEN, muB=MU_B_DEFAULT,
                      p=P_DEPH, steps=36,
                      savepath="entropy_R_page_like.png")

    # Print absolute paths so you can open them directly
    for f in ["MI_RH_RS_vs_muB.png", "MI_RH_RS_vs_p_for_muB.png", "entropy_R_page_like.png"]:
        print("Saved:", os.path.abspath(f))
