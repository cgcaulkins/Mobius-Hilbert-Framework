
Möbius Decoder Sim

mobius_decoder_sim.py

Dr. Chris Caulkins
c.caulkins@FreeRangeThinker.org

Overview

This Python script is a numerical simulation designed to test the decodability of information within the Möbius-Hilbert framework. It models a tripartite quantum system consisting of a Source (S), a Hidden domain (H), and Radiation (R).

The core purpose of this simulation is to answer the question: "Where does the information live?"

It does this by encoding a secret message (a quantum phase, θ) into the system and then calculating an observer's theoretical ability to decode that message. The observer's success is measured by the trace distance, which quantifies how distinguishable two different messages are.

    A trace distance of 0 means the messages are indistinguishable (no information is learnable).

    A trace distance > 0 means the messages are distinguishable (information is learnable).

The simulation compares the decodability for three different observers:

    An observer who can only see the Radiation (R).

    An observer who can see the Radiation and the Source (R ∪ S).

    An observer who can see the Radiation and the Hidden domain (R ∪ H).

Core Concepts Modeled

    step_density(mu, theta, p): This function creates the central quantum state for a single emission. It models the core hypothesis: the information is in a superposition of being in the Source-Radiation subsystem (with probability 1−μ) and the Hidden-Radiation subsystem (with probability μ).

    dephase_on_R(rho, p): This applies an observation/noise channel to the Radiation qubit, simulating the effect of environmental decoherence (p) on the escaping Hawking particle.

    run_experiment(...): This function calculates how the trace distance (decodability) grows as an observer collects more and more emissions (k).

    run_and_plot(): This main function runs the experiment for three hard-coded scenarios to generate the key figures for the paper.

Requirements

The script requires the following Python libraries:

    numpy (for quantum state linear algebra)

    matplotlib (for generating the plots)

    pathlib (for handling file paths, included in standard Python 3)

You can install the dependencies using pip:

pip install numpy matplotlib

How to Run

You can run the script directly from your terminal. It will automatically execute the run_and_plot() function, which runs all three scenarios and saves the resulting plots in the same directory.

python3 mobius_decoder_sim.py

Outputs

This script will generate and save six PNG files corresponding to the three scenarios:

    Scenario: μ=0.0,p=0.0 (No migration, no noise)

        radiation_trace_mu0.0_p0.0.png (Shows trace distance for R only)

        paired_trace_mu0.0_p0.0.png (Compares R ∪ S vs. R ∪ H)

    Scenario: μ=0.6,p=0.0 (With migration, no noise)

        radiation_trace_mu0.6_p0.0.png (Shows trace distance for R only)

        paired_trace_mu0.6_p0.0.png (Compares R ∪ S vs. R ∪ H)

    Scenario: μ=0.6,p=0.5 (With migration, with noise)

        radiation_trace_mu0.6_p0.5.png (Shows trace distance for R only)

        paired_trace_mu0.6_p0.5.png (Compares R ∪ S vs. R ∪ H)