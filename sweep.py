import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------
# Config
# -------------------------------
DATA_CSV = "bh_data.csv"   # change if your file is named differently

# Helper: ensure directory exists
def _ensure_dir(path):
    d = os.path.dirname(path)
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)

# -------------------------------
# Data loading with graceful fallback
# -------------------------------
def load_data(csv_path):
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        # Normalize potential tidy format
        if "metric" in df.columns and "value" in df.columns:
            # Pivot long -> wide by metric
            keys = [c for c in ["mu","p","RH","RS","SH"] if c in df.columns]
            df_wide = df.pivot_table(index=keys, columns="metric", values="value", aggfunc="mean").reset_index()
            df_wide.columns.name = None
            return df_wide
        return df
    # Synthetic fallback
    rng = np.random.default_rng(42)
    mus = np.linspace(0.2, 1.0, 21)
    ps  = np.linspace(0.1, 0.9, 17)
    RHs = np.linspace(0.1, 1.0, 20)
    RSs = np.linspace(0.1, 1.0, 20)
    SHs = np.linspace(0.1, 1.0, 25)

    # Simple toy functions for demonstration (replace with your physics/estimates)
    def f_decoder(mu, p):  # peaks around mid p & higher mu
        return 0.5 + 0.45*np.tanh(2.0*(mu-0.5))*np.exp(-((p-0.5)**2)/0.05)

    def f_MI_mu(mu, p):    # grows with mu, modulated by p
        return 0.1 + 0.9*(mu**1.2)*(0.6+0.4*np.cos(2*np.pi*p))

    def f_MI_RH(RH, p):    # increases with RH and p
        return 0.05 + 0.95*(RH**0.7)*(0.4+0.6*p)

    def f_MI_RS(RS, p):    # different slope vs RS
        return 0.03 + 0.85*(RS**0.9)*(0.5+0.5*p)

    def f_MI_SH(SH, p):    # saturating in SH; depends on p
        return 0.02 + (1.0 - np.exp(-3*SH))*(0.3+0.7*p)

    rows = []
    for mu in mus:
        for p in ps:
            rows.append({
                "mu": mu, "p": p,
                "decoder_acc": f_decoder(mu,p),
                "MI_mu": f_MI_mu(mu,p)
            })
    for p in ps:
        for RH, RS in zip(RHs, RSs):
            rows.append({"p": p, "RH": RH, "MI_RH": f_MI_RH(RH,p)})
            rows.append({"p": p, "RS": RS, "MI_RS": f_MI_RS(RS,p)})
        for SH in SHs:
            rows.append({"p": p, "SH": SH, "MI_SH": f_MI_SH(SH,p)})
    df = pd.DataFrame(rows).drop_duplicates()
    return df

df = load_data(DATA_CSV)

# -------------------------------
# 1) decoder_accuracy_mu0.6_p0.5.png
# -------------------------------
def plot_decoder_accuracy_fixed_mu_p(df, mu_fixed=0.6, p_fixed=0.5, outpath="decoder_accuracy_mu0.6_p0.5.png"):
    # Try exact match. If not, interpolate over a mesh created from df.
    subset = df.copy()
    if "decoder_acc" not in subset.columns:
        print("[warn] decoder_acc not found; skipping this figure.")
        return
    # If we have exact rows for mu,p:
    exact = subset[(subset.get("mu").notna()) & (subset.get("p").notna())]
    exact = exact[(np.isclose(exact["mu"], mu_fixed)) & (np.isclose(exact["p"], p_fixed))]
    _ensure_dir(outpath)
    if not exact.empty:
        val = exact["decoder_acc"].iloc[0]
        plt.figure()
        plt.title(f"Decoder Accuracy at μ={mu_fixed}, p={p_fixed}")
        plt.bar(["decoder_acc"], [val])
        plt.ylim(0, 1.0)
        plt.ylabel("Accuracy")
        plt.tight_layout()
        plt.savefig(outpath, dpi=200)
        plt.close()
        return

    # Otherwise, build a 2D map and show a local slice around (mu_fixed,p_fixed)
    if all(c in subset.columns for c in ["mu","p","decoder_acc"]):
        # nearest neighbors around target point
        subset["dist2"] = (subset["mu"]-mu_fixed)**2 + (subset["p"]-p_fixed)**2
        nearest = subset.nsmallest(200, "dist2")
        plt.figure()
        plt.title(f"Decoder Accuracy near μ={mu_fixed}, p={p_fixed}")
        plt.scatter(nearest["mu"], nearest["p"], s=20, alpha=0.7, label="samples")
        plt.xlabel("μ")
        plt.ylabel("p")
        sc = plt.scatter([mu_fixed],[p_fixed], s=80, marker="x", label="target")
        # annotate estimated local value
        est = nearest.iloc[0]["decoder_acc"]
        plt.text(mu_fixed, p_fixed, f" ~{est:.3f}", ha="left", va="bottom")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outpath, dpi=200)
        plt.close()
    else:
        print("[warn] Not enough columns to estimate decoder_acc around (mu,p).")

plot_decoder_accuracy_fixed_mu_p(df, 0.6, 0.5, "decoder_accuracy_mu0.6_p0.5.png")

# -------------------------------
# 2) MI_mu_sweep.png  (MI vs mu, for a few representative p values)
# -------------------------------
def plot_MI_mu_sweep(df, p_values=(0.2,0.5,0.8), outpath="MI_mu_sweep.png"):
    cols_ok = {"mu","p","MI_mu"}.issubset(df.columns)
    if not cols_ok:
        print("[warn] Columns for MI_mu_sweep not found; skipping.")
        return
    plt.figure()
    plt.title("Mutual Information (MI) vs μ (selected p)")
    for pv in p_values:
        sl = df[(df["p"].notna()) & (np.isclose(df["p"], pv))]
        sl = sl[sl["mu"].notna()].sort_values("mu")
        if sl.empty: 
            continue
        plt.plot(sl["mu"], sl["MI_mu"], label=f"p={pv}")
    plt.xlabel("μ")
    plt.ylabel("MI(μ; •)")
    plt.legend()
    plt.tight_layout()
    _ensure_dir(outpath)
    plt.savefig(outpath, dpi=200)
    plt.close()

plot_MI_mu_sweep(df, (0.2,0.5,0.8), "MI_mu_sweep.png")

# -------------------------------
# 3) MI_RH_vs_RS_sweep_p.png (MI_RH and MI_RS vs horizon radii across p)
# -------------------------------
def plot_MI_RH_RS_sweep_p(df, p_values=(0.2,0.5,0.8), outpath="MI_RH_vs_RS_sweep_p.png"):
    need_cols_RH = {"RH","p","MI_RH"}.issubset(df.columns)
    need_cols_RS = {"RS","p","MI_RS"}.issubset(df.columns)
    if not (need_cols_RH or need_cols_RS):
        print("[warn] Columns for MI_RH/MI_RS not found; skipping.")
        return
    plt.figure()
    plt.title("MI vs Horizon Radii (RH, RS) across p")
    # RH
    if need_cols_RH:
        for pv in p_values:
            sl = df[(df["p"].notna()) & (np.isclose(df["p"], pv)) & df["RH"].notna()]
            sl = sl.sort_values("RH")
            if not sl.empty:
                plt.plot(sl["RH"], sl["MI_RH"], linestyle="-", label=f"MI_RH, p={pv}")
    # RS
    if need_cols_RS:
        for pv in p_values:
            sl = df[(df["p"].notna()) & (np.isclose(df["p"], pv)) & df["RS"].notna()]
            sl = sl.sort_values("RS")
            if not sl.empty:
                plt.plot(sl["RS"], sl["MI_RS"], linestyle="--", label=f"MI_RS, p={pv}")
    plt.xlabel("Horizon radius (R)")
    plt.ylabel("MI")
    plt.legend()
    plt.tight_layout()
    _ensure_dir(outpath)
    plt.savefig(outpath, dpi=200)
    plt.close()

plot_MI_RH_RS_sweep_p(df, (0.2,0.5,0.8), "MI_RH_vs_RS_sweep_p.png")

# -------------------------------
# 4) MI_SH_sweep_p.png (MI_SH vs SH across p)
# -------------------------------
def plot_MI_SH_sweep_p(df, p_values=(0.2,0.5,0.8), outpath="MI_SH_sweep_p.png"):
    if not {"SH","p","MI_SH"}.issubset(df.columns):
        print("[warn] Columns for MI_SH not found; skipping.")
        return
    plt.figure()
    plt.title("MI vs Surface Gravity/Entropy Proxy (SH) across p")
    for pv in p_values:
        sl = df[(df["p"].notna()) & (np.isclose(df["p"], pv)) & df["SH"].notna()]
        sl = sl.sort_values("SH")
        if not sl.empty:
            plt.plot(sl["SH"], sl["MI_SH"], label=f"p={pv}")
    plt.xlabel("SH")
    plt.ylabel("MI_SH")
    plt.legend()
    plt.tight_layout()
    _ensure_dir(outpath)
    plt.savefig(outpath, dpi=200)
    plt.close()

plot_MI_SH_sweep_p(df, (0.2,0.5,0.8), "MI_SH_sweep_p.png")

print("Done. Wrote figures if data available (or synthetic fallback used).")
