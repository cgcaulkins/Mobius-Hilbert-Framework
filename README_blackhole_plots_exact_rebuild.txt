
Black Hole Plots — Exact Rebuild Instructions
============================================

Goal
----
Recreate the four figures with *identical* filenames and shapes as your last program:

- decoder_accuracy_mu0.6_p0.5.png
- MI_mu_sweep.png
- MI_RH_vs_RS_sweep_p.png
- MI_SH_sweep_p.png

What you need to do
-------------------
1) Open `blackhole_plots_exact_rebuild.py` and paste your ORIGINAL equations into the
   section marked:
       ==== PASTE YOUR ORIGINAL EQUATIONS HERE ====

   Replace the NotImplementedError blocks in:
   - decoder_accuracy(mu, p, row=None)
   - MI_mu(mu, p, row=None)
   - MI_RH(RH, p, row=None)
   - MI_RS(RS, p, row=None)
   - MI_SH(SH, p, row=None)

   Keep the function signatures. If your prior code referenced additional fields,
   you can access them via the `row` argument (a pandas namedtuple of the current row).

2) Provide your data CSV as `bh_data.csv` (or set BH_DATA_CSV env var). The script accepts:
   a) Wide format (one row per parameter combo):
      mu, p, RH, RS, SH, <any extra fields your equations require>
      Optionally you can include precomputed columns if your equations used them.

   b) Tidy/long format:
      metric, mu, p, RH, RS, SH, value
      (The script pivots this to wide automatically.)

   A small template CSV is included as `bh_data_template.csv` for column guidance.

3) Run:
   python blackhole_plots_exact_rebuild.py
   (Optionally set output directory: BH_OUT_DIR=/path/to/out)

Matplotlib parity
-----------------
For bit-for-bit identical images, try to match your previous matplotlib version
and rcParams. If your old script set fonts, linewidths, or styles, paste those
lines near the top of this script to ensure perfect visual parity.

Troubleshooting
---------------
- KeyError: 'SH' (or others): ensure your CSV has that column or adapt the code to your schema.
- Wrong numeric values: verify the pasted equations and that the CSV matches the one used last time.
- Slight visual differences: harmonize matplotlib versions and any style/rc settings you used before.
