
Parameter Sweep for Black Hole Simulations 

sweep.py

Dr. Chris Caulkins
c.caulkins@FreeRangeThinker.org

Overview

This Python script, sweep.py, is a utility for running parameter sweeps of the quantum simulations related to the Möbius-Hilbert framework. It is designed to be used with a primary simulation script (like mobius_hypermag_model.py or a similar script) to test the model's behavior across a range of different physical parameters.

The script automates the process of running the simulation multiple times with different values for key parameters and collecting the results. This is essential for understanding how the model's predictions (like mutual information and entropy) depend on factors such as the migration ratio (mu) and dephasing strength (p).

Purpose

The main goal of this script is to generate the data needed to create plots that show trends and relationships between the model's parameters and its outputs. For example, it can be used to generate the data for a plot of "Mutual Information vs. Migration Rate (mu)" by running the simulation for many different mu values.

How It Works

The script is structured to:

    Define ranges for the parameters to be tested (e.g., a list of mu values and p values).

    Loop through these parameter combinations.

    For each combination, it calls the main simulation function (e.g., simulate()) from another script.

    It collects the output data from each simulation run.

    Finally, it uses this collected data to generate plots that visualize the results of the   parameter sweep.

How to Use

This script is not intended to be run as a standalone file. It is a tool to be used in conjunction with one of the primary simulation scripts. To use it, you would typically:

    Import the sweep.py script into your main simulation file or a Jupyter Notebook.

    Call the functions defined in sweep.py to run your desired parameter sweeps.

    The script will then generate and save the corresponding plots (like MI_mu_sweep.png).

Note: The specific implementation details (like which simulation function is called and what parameters are swept) are defined within the script's code.

Requirements

This script has the same dependencies as the main simulation scripts:

    numpy

    matplotlib

    pandas

You can install these dependencies using pip:
pip install numpy pandas matplotlib

How to Run

This script is designed to be run after you have run your primary simulations and collected their output into a single CSV file.

Step 1: Prepare Your Data

This script requires a data file named bh_data.csv to be in the same directory. You must create this file manually by running your other simulation scripts (like mobius_hypermag_model.py, mobius_decoder_sim.py, etc.) and compiling their numerical outputs into a single CSV.

The script expects bh_data.csv to have columns like mu, p, decoder_acc, MI_mu, RH, MI_RH, etc.

tep 2: Run the Script

Once bh_data.csv is in the same folder, run the script from your terminal:
python3 sweep.py

he script will load bh_data.csv, generate the plots, and save them as PNG files in the same directory.

Note on Synthetic Data

If you run the script without a bh_data.csv file, it will not fail. It will automatically generate a set of synthetic, "toy" data to plot instead. The resulting plots will look reasonable but will not be based on your actual simulation results.

Outputs

This script will generate and save the following four PNG files:

    decoder_accuracy_mu0.6_p0.5.png: A plot visualizing decoder accuracy, typically for a fixed set of parameters.

    MI_mu_sweep.png: A plot showing how Mutual Information changes as the migration rate (μ) is swept, for a few different dephasing (p) values.

    MI_RH_vs_RS_sweep_p.png: A comparison plot of I(R:H) and I(R:S) against their respective parameters.

    MI_SH_sweep_p.png: A plot showing the mutual information between the Source (S) and Hidden (H) domains.

