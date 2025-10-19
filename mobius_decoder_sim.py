# mobius_decoder_sim.py
import numpy as np
import matplotlib.pyplot as plt
from functools import reduce
from pathlib import Path

# ---------- helpers ----------
def kron_all(mats):
    return reduce(lambda A,B: np.kron(A,B), mats)

def to_dm(psi):
    return psi @ psi.conj().T

def dephase_on_R(rho, p):
    """Apply Schrödinger-style observation (dephasing) to the R qubit of [S,H,R]."""
    K0 = np.sqrt(1-p) * np.eye(2, dtype=complex)
    K1 = np.sqrt(p)   * np.array([[1,0],[0,0]], dtype=complex)
    K2 = np.sqrt(p)   * np.array([[0,0],[0,1]], dtype=complex)
    I = np.eye(2, dtype=complex)
    out = np.zeros_like(rho, dtype=complex)
    for K in (K0,K1,K2):
        U = np.kron(np.kron(I, I), K)  # act on R only
        out += U @ rho @ U.conj().T
    return out

# ---------- single emission (isometry slice) ----------
def step_density(mu, theta, p):
    """
    Build one emission state over 3 qubits ordered [S, H, R]:
      sqrt(1-μ) |Φ_θ>_{S,R}|0>_H  +  sqrt(μ) |0>_S |Φ_θ>_{H,R}
    with |Φ_θ> = (|00> + e^{iθ}|11>)/√2
    Then apply dephasing of strength p on R (observation).
    """
    a = np.sqrt(1-mu)
    b = np.sqrt(mu)
    phase = np.exp(1j*theta)
    v = np.zeros(8, dtype=complex)  # basis |s h r>
    # a * |Φ_θ>_{S,R}|0>_H  contributes |000> and |101>
    v[0] += a/np.sqrt(2)          # |0 0 0>
    v[5] += a*phase/np.sqrt(2)    # |1 0 1>
    # b * |0>_S |Φ_θ>_{H,R}  contributes |000> and |011>
    v[0] += b/np.sqrt(2)          # |0 0 0> (coherent overlap)
    v[3] += b*phase/np.sqrt(2)    # |0 1 1>
    v /= np.linalg.norm(v)
    rho = to_dm(v.reshape(8,1))
    if p > 0:
        rho = dephase_on_R(rho, p)
    return rho  # density matrix over [S,H,R]

def reduce_single_step(rho_SHR):
    """Return reduced ρ_R (2x2), ρ_RS (4x4), ρ_RH (4x4) for one [S,H,R] step."""
    rho6 = rho_SHR.reshape(2,2,2, 2,2,2)  # (S,H,R ; S',H',R')

    # ρ_R[r,r'] = sum_{s,h} ρ[s,h,r ; s,h,r']
    rho_R = np.zeros((2,2), dtype=complex)
    for r in range(2):
        for rp in range(2):
            for S in range(2):
                for H in range(2):
                    rho_R[r,rp] += rho6[S,H,r, S,H,rp]

    # ρ_RS[(s,r),(s',r')] = sum_h ρ[s,h,r ; s',h,r']
    rho_RS = np.zeros((4,4), dtype=complex)
    for S in range(2):
        for R in range(2):
            for Sp in range(2):
                for Rp in range(2):
                    val = 0.0+0j
                    for H in range(2):
                        val += rho6[S,H,R, Sp,H,Rp]
                    i = S*2 + R
                    j = Sp*2 + Rp
                    rho_RS[i,j] = val

    # ρ_RH[(h,r),(h',r')] = sum_s ρ[s,h,r ; s,h',r']
    rho_RH = np.zeros((4,4), dtype=complex)
    for H in range(2):
        for R in range(2):
            for Hp in range(2):
                for Rp in range(2):
                    val = 0.0+0j
                    for S in range(2):
                        val += rho6[S,H,R, S,Hp,Rp]
                    i = H*2 + R
                    j = Hp*2 + Rp
                    rho_RH[i,j] = val

    return rho_R, rho_RS, rho_RH

def trace_distance(rho, sigma):
    """0.5 * ||rho - sigma||_1; for Hermitian diff, sum absolute eigenvalues."""
    diff = (rho - sigma + (rho - sigma).conj().T)/2.0
    eigs = np.linalg.eigvalsh(diff)
    return 0.5*float(np.sum(np.abs(np.real(eigs))))

# ---------- experiment driver ----------
def run_experiment(mu, p, n_R=6, n_pairs=4, theta1=np.pi/2):
    """
    Compare two histories (θ=0 vs θ=θ1). Build per-step reduced states,
    then tensor-prefix them to see how distinguishability grows with emissions.
    """
    # Precompute one-step reductions for θ=0 and θ=θ1
    R0_list=[]; RS0_list=[]; RH0_list=[]
    R1_list=[]; RS1_list=[]; RH1_list=[]
    for _ in range(max(n_R, n_pairs)):
        rho0 = step_density(mu, 0.0, p)
        rho1 = step_density(mu, theta1, p)
        r0, rs0, rh0 = reduce_single_step(rho0)
        r1, rs1, rh1 = reduce_single_step(rho1)
        R0_list.append(r0); RS0_list.append(rs0); RH0_list.append(rh0)
        R1_list.append(r1); RS1_list.append(rs1); RH1_list.append(rh1)

    # Grow prefixes: R only
    kR=[]; TD_R=[]; rhoR0=None; rhoR1=None
    for k in range(1, n_R+1):
        rhoR0 = R0_list[k-1] if rhoR0 is None else np.kron(rhoR0, R0_list[k-1])
        rhoR1 = R1_list[k-1] if rhoR1 is None else np.kron(rhoR1, R1_list[k-1])
        TD_R.append(trace_distance(rhoR0, rhoR1)); kR.append(k)

    # Grow prefixes: RS and RH (paired views)
    kP=[]; TD_RS=[]; TD_RH=[]; rhoRS0=None; rhoRS1=None; rhoRH0=None; rhoRH1=None
    for k in range(1, n_pairs+1):
        rhoRS0 = RS0_list[k-1] if rhoRS0 is None else np.kron(rhoRS0, RS0_list[k-1])
        rhoRS1 = RS1_list[k-1] if rhoRS1 is None else np.kron(rhoRS1, RS1_list[k-1])
        rhoRH0 = RH0_list[k-1] if rhoRH0 is None else np.kron(rhoRH0, RH0_list[k-1])
        rhoRH1 = RH1_list[k-1] if rhoRH1 is None else np.kron(rhoRH1, RH1_list[k-1])
        TD_RS.append(trace_distance(rhoRS0, rhoRS1))
        TD_RH.append(trace_distance(rhoRH0, rhoRH1))
        kP.append(k)

    return (np.array(kR), np.array(TD_R), np.array(kP), np.array(TD_RS), np.array(TD_RH))

# ---------- plotting ----------
def run_and_plot():
    outdir = Path(".")
    scenarios = [
        {"mu":0.0, "p":0.0, "label":"μ=0.0, p=0.0"},
        {"mu":0.6, "p":0.0, "label":"μ=0.6, p=0.0"},
        {"mu":0.6, "p":0.5, "label":"μ=0.6, p=0.5"},
    ]
    figs = []
    for sc in scenarios:
        kR, TD_R, kP, TD_RS, TD_RH = run_experiment(mu=sc["mu"], p=sc["p"], n_R=6, n_pairs=4, theta1=np.pi/2)

        # Radiation only
        fig1 = plt.figure(figsize=(7,5))
        plt.plot(kR, TD_R, marker='o')
        plt.xlabel("Emissions k (R only)")
        plt.ylabel("Trace distance (θ=0 vs θ=π/2)")
        plt.title(f"Learnable signal in radiation — {sc['label']}")
        plt.tight_layout()
        fig1.savefig(outdir / f"radiation_trace_{sc['label'].replace(' ','_').replace('μ','mu')}.png", dpi=150)
        figs.append(fig1)

        # Where the message lives: R∪S vs R∪H
        fig2 = plt.figure(figsize=(7,5))
        plt.plot(kP, TD_RS, marker='o', label="R ∪ S")
        plt.plot(kP, TD_RH, marker='s', label="R ∪ H")
        plt.xlabel("Emissions k (paired views)")
        plt.ylabel("Trace distance (θ=0 vs θ=π/2)")
        plt.title(f"Where the message lives — {sc['label']}")
        plt.legend()
        plt.tight_layout()
        fig2.savefig(outdir / f"paired_trace_{sc['label'].replace(' ','_').replace('μ','mu')}.png", dpi=150)
        figs.append(fig2)
    return figs

if __name__ == "__main__":
    run_and_plot()
    # Show windows so they don't vanish when the script exits
    plt.show()
Xx