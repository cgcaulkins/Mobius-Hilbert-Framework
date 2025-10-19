import numpy as np
import matplotlib.pyplot as plt

zero = np.array([[1.0],[0.0]], dtype=complex)
one  = np.array([[0.0],[1.0]], dtype=complex)

def build_step_density(mu):
    a = np.sqrt(1-mu)
    b = np.sqrt(mu)
    N = 1/np.sqrt(1 + a*b)
    v = np.zeros(8, dtype=complex)  # (S,H,R) indexing as binary 3-bit number
    v[0] += N*(a/np.sqrt(2) + b/np.sqrt(2))  # |000>
    v[5] += N*(a/np.sqrt(2))                 # |101>
    v[3] += N*(b/np.sqrt(2))                 # |011>
    psi = v.reshape(8,1)
    rho = psi @ psi.conj().T
    return rho

def dephasing_on_R(rho, p):
    # Kraus on R only
    K0 = np.sqrt(1-p) * np.eye(2, dtype=complex)
    K1 = np.sqrt(p) * np.array([[1,0],[0,0]], dtype=complex)
    K2 = np.sqrt(p) * np.array([[0,0],[0,1]], dtype=complex)
    I = np.eye(2, dtype=complex)
    out = np.zeros_like(rho, dtype=complex)
    for K in (K0, K1, K2):
        U = np.kron(np.kron(I, I), K)  # S,H,R
        out += U @ rho @ U.conj().T
    return out

def reduce_R_and_SH(rho):
    # reshape to (S,H,R,S',H',R')
    rho6 = rho.reshape(2,2,2, 2,2,2)
    # rho_R[r,r'] = sum_{s,h} rho6[s,h,r, s,h,r']
    rho_R = np.zeros((2,2), dtype=complex)
    for r in range(2):
        for rp in range(2):
            s = 0
            for S in range(2):
                for H in range(2):
                    rho_R[r,rp] += rho6[S,H,r, S,H,rp]
    # rho_SH[s,h ; s',h'] = sum_{r} rho6[s,h,r, s',h',r]
    rho_SH = np.zeros((4,4), dtype=complex)
    for S in range(2):
        for H in range(2):
            for Sp in range(2):
                for Hp in range(2):
                    val = 0.0+0j
                    for r in range(2):
                        val += rho6[S,H,r, Sp,Hp,r]
                    idx = S*2 + H
                    jdx = Sp*2 + Hp
                    rho_SH[idx,jdx] = val
    return rho_R, rho_SH

def vN_entropy(rho, tol=1e-12):
    rho = (rho + rho.conj().T)/2.0
    vals = np.linalg.eigvalsh(rho)
    vals = np.real(vals)
    vals = vals[vals>tol]
    return float(-np.sum(vals*np.log2(vals)))

def step_entropies(mu, p):
    rho = build_step_density(mu)
    if p>0:
        rho = dephasing_on_R(rho, p)
    rho_R, rho_SH = reduce_R_and_SH(rho)
    SR = vN_entropy(rho_R)
    SSH = vN_entropy(rho_SH)
    STOT = vN_entropy(rho)
    I = SR + SSH - STOT
    return SR, I

def simulate(Nsteps=32, mu=0.3, p=0.0):
    S_cum = []
    I_cum = []
    ssum = 0.0
    isum = 0.0
    for k in range(1, Nsteps+1):
        SR, I = step_entropies(mu, p)
        ssum += SR
        isum += I
        S_cum.append(ssum)
        I_cum.append(isum)
    return np.arange(1, Nsteps+1), np.array(S_cum), np.array(I_cum)

# Run and plot
Nsteps = 32
mus = [0.0, 0.3, 0.6, 0.9]
p_values = [0.0, 0.5]

for p in p_values:
    k = None
    plt.figure(figsize=(7,5))
    for mu in mus:
        k, Srad, Icum = simulate(Nsteps=Nsteps, mu=mu, p=p)
        plt.plot(k, Srad, label=f"μ={mu:.1f}")
    plt.xlabel("Emission count k")
    plt.ylabel("Cumulative radiation entropy S_R(k) [qubits]")
    plt.title(f"Radiation entropy with observation strength p={p}")
    plt.legend()
    plt.tight_layout()

for p in p_values:
    plt.figure(figsize=(7,5))
    for mu in mus:
        k, Srad, Icum = simulate(Nsteps=Nsteps, mu=mu, p=p)
        plt.plot(k, Icum, label=f"μ={mu:.1f}")
    plt.xlabel("Emission count k")
    plt.ylabel("Cumulative I(R : S∪H) [bits]")
    plt.title(f"Radiation–(Singularity+Hidden) mutual information, p={p}")
    plt.legend()
    plt.tight_layout()
