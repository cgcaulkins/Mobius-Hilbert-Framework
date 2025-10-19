#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Black Hole Research: Exact Plot Rebuild Scaffold
------------------------------------------------
This script is structured to reproduce the four figures from your prior program,
with the same filenames and plot organization. Paste your *original equations*
into the clearly marked TODO section below to make the outputs bit-for-bit
identical (assuming the same data source and matplotlib version).

Expected outputs (written to the current working directory):
- decoder_accuracy_mu0.6_p0.5.png
- MI_mu_sweep.png
- MI_RH_vs_RS_sweep_p.png
- MI_SH_sweep_p.png

Data handling:
- Works with either "wide" or "tidy/long" CSV formats.
- See README for the column schema.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------
# Config
# -------------------------------
DATA_CSV = os.environ.get("BH_DATA_CSV", "bh_data.csv")
OUT_DIR = os.environ.get("BH_OUT_DIR", ".")

# -------------------------------
# ==== PASTE YOUR ORIGINAL EQUATIONS HERE ====
# Replace the NotImplementedError blocks with your prior program's logic.
# Keep the signatures intact so the plotting code below "just works".
# If your original code used different variable names, map them here.
# -------------------------------

def decoder_accuracy(mu: float, p: float, row=None) -> float:
    """
    Return decoder accuracy for the given mu, p.
    If you computed from per-row features, use `row` (a pandas Series).
    """
    raise NotImplementedError(
        "TODO: Paste your original decoder accuracy equation here."
    )

def MI_mu(mu: float, p: float, row=None) -> float:
    """
    Mutual information as a function of mu and (optionally) p.
    """
    raise NotImplementedError(
        "TODO: Paste your original MI(mu, p) equation here."
    )

def MI_RH(RH: float, p: float, row=None) -> float:
    """
    Mutual information vs RH for a given p.
    """
    raise NotImplementedError(
        "TODO: Paste your original MI_RH(RH, p) equation here."
    )

def MI_RS(RS: float, p: float, row=None) -> float:
    """
    Mutual information vs RS for a given p.
    """
    raise NotImplementedError(
        "TODO: Paste your original MI_RS(RS, p) equation here."
    )

def MI_SH(SH: float, p: float, row=None) -> float:
    """
    Mutual information vs SH for a given p.
    """
    raise NotImplementedError(
        "TODO: Paste your original MI_SH(SH, p) equation here."
    )

# -------------------------------
# Data loading & normalization
# -------------------------------

def load_data(csv_path: str) -> pd.DataFrame:
    if not os.path.exists(csv_path):
        raise FileNotFoundError(
            f"Data file '{csv_path}' not found. Set BH_DATA_CSV or place CSV next to the script."
        )
    df = pd.read_csv(csv_path)
    # Normalize potential tidy format
    if "metric" in df.columns and "value" in df.columns:
        keys = [c for c in ["mu","p","RH","RS","SH"] if c in df.columns]
        df_wide = df.pivot_table(index=keys, columns="metric", values="value", aggfunc="mean").reset_index()
        df_wide.columns.name = None
        return df_wide
    return df

def _ensure_dir(path: str):
    d = os.path.dirname(path)
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)

# -------------------------------
# Plot 1: decoder_accuracy_mu0.6_p0.5.png
# -------------------------------
def plot_decoder_accuracy_fixed_mu_p(df: pd.DataFrame, mu_fixed=0.6, p_fixed=0.5, fname="decoder_accuracy_mu0.6_p0.5.png"):
    outpath = os.path.join(OUT_DIR, fname)
    _ensure_dir(outpath)

    # Try to find a matching row for exact evaluation
    has_mu = "mu" in df.columns
    has_p = "p" in df.columns
    val = None
    if has_mu and has_p:
        exact = df[np.isclose(df["mu"], mu_fixed) & np.isclose(df["p"], p_fixed)]
        if not exact.empty:
            # If original code used row-dependent fields, pass row
            val = decoder_accuracy(mu_fixed, p_fixed, row=exact.iloc[0])
        else:
            # Fallback: nearest row
            tmp = df.copy()
            tmp["dist2"] = (tmp.get("mu", mu_fixed) - mu_fixed)**2 + (tmp.get("p", p_fixed) - p_fixed)**2
            nearest = tmp.nsmallest(1, "dist2")
            if not nearest.empty:
                val = decoder_accuracy(nearest.iloc[0].get("mu", mu_fixed), nearest.iloc[0].get("p", p_fixed), row=nearest.iloc[0])
    if val is None:
        # Last resort: direct function without row
        val = decoder_accuracy(mu_fixed, p_fixed, row=None)

    plt.figure()
    plt.title(f"Decoder Accuracy at μ={mu_fixed}, p={p_fixed}")
    plt.bar(["decoder_acc"], [val])
    plt.ylim(0, 1.0)
    plt.ylabel("Accuracy")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

# -------------------------------
# Plot 2: MI_mu_sweep.png
# -------------------------------
def plot_MI_mu_sweep(df: pd.DataFrame, p_values=(0.2,0.5,0.8), fname="MI_mu_sweep.png"):
    outpath = os.path.join(OUT_DIR, fname)
    _ensure_dir(outpath)

    plt.figure()
    plt.title("Mutual Information (MI) vs μ (selected p)")
    for pv in p_values:
        sl = df.copy()
        if "p" in sl.columns:
            sl = sl[np.isclose(sl["p"], pv)]
        if "mu" in sl.columns:
            sl = sl.dropna(subset=["mu"]).sort_values("mu")
            if not sl.empty:
                y = [MI_mu(mu=float(m), p=pv, row=row) for m, row in zip(sl["mu"].values, sl.itertuples(index=False))]
                plt.plot(sl["mu"], y, label=f"p={pv}")
    plt.xlabel("μ")
    plt.ylabel("MI(μ; •)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

# -------------------------------
# Plot 3: MI_RH_vs_RS_sweep_p.png
# -------------------------------
def plot_MI_RH_RS_sweep_p(df: pd.DataFrame, p_values=(0.2,0.5,0.8), fname="MI_RH_vs_RS_sweep_p.png"):
    outpath = os.path.join(OUT_DIR, fname)
    _ensure_dir(outpath)

    plt.figure()
    plt.title("MI vs Horizon Radii (RH, RS) across p")

    if "RH" in df.columns:
        for pv in p_values:
            sl = df.copy()
            if "p" in sl.columns:
                sl = sl[np.isclose(sl["p"], pv)]
            sl = sl.dropna(subset=["RH"]).sort_values("RH")
            if not sl.empty:
                y = [MI_RH(RH=float(rh), p=pv, row=row) for rh, row in zip(sl["RH"].values, sl.itertuples(index=False))]
                plt.plot(sl["RH"], y, linestyle="-", label=f"MI_RH, p={pv}")

    if "RS" in df.columns:
        for pv in p_values:
            sl = df.copy()
            if "p" in sl.columns:
                sl = sl[np.isclose(sl["p"], pv)]
            sl = sl.dropna(subset=["RS"]).sort_values("RS")
            if not sl.empty:
                y = [MI_RS(RS=float(rs), p=pv, row=row) for rs, row in zip(sl["RS"].values, sl.itertuples(index=False))]
                plt.plot(sl["RS"], y, linestyle="--", label=f"MI_RS, p={pv}")

    plt.xlabel("Horizon radius (R)")
    plt.ylabel("MI")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

# -------------------------------
# Plot 4: MI_SH_sweep_p.png
# -------------------------------
def plot_MI_SH_sweep_p(df: pd.DataFrame, p_values=(0.2,0.5,0.8), fname="MI_SH_sweep_p.png"):
    outpath = os.path.join(OUT_DIR, fname)
    _ensure_dir(outpath)

    if "SH" not in df.columns:
        raise KeyError("Column 'SH' not found in data.")

    plt.figure()
    plt.title("MI vs Surface Gravity/Entropy Proxy (SH) across p")
    for pv in p_values:
        sl = df.copy()
        if "p" in sl.columns:
            sl = sl[np.isclose(sl["p"], pv)]
        sl = sl.dropna(subset=["SH"]).sort_values("SH")
        if not sl.empty:
            y = [MI_SH(SH=float(sh), p=pv, row=row) for sh, row in zip(sl["SH"].values, sl.itertuples(index=False))]
            plt.plot(sl["SH"], y, label=f"p={pv}")
    plt.xlabel("SH")
    plt.ylabel("MI_SH")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

# -------------------------------
# Main
# -------------------------------
def main():
    df = load_data(DATA_CSV)
    plot_decoder_accuracy_fixed_mu_p(df, 0.6, 0.5, "decoder_accuracy_mu0.6_p0.5.png")
    plot_MI_mu_sweep(df, (0.2,0.5,0.8), "MI_mu_sweep.png")
    plot_MI_RH_RS_sweep_p(df, (0.2,0.5,0.8), "MI_RH_vs_RS_sweep_p.png")
    plot_MI_SH_sweep_p(df, (0.2,0.5,0.8), "MI_SH_sweep_p.png")
    print("Done. Figures written to:", os.path.abspath(OUT_DIR))

if __name__ == "__main__":
    main()
