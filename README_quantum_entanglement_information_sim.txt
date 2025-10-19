
Quantum Entanglement Information Sim   

quantum_entanglement_information_sim.py

Dr. Chris Caulkins
c.caulkins@FreeRangeThinker.org

Overview

This Python script is a "toy model" simulation designed to illustrate the cumulative information and entropy dynamics of the Möbius-Hilbert framework. It models a 3-qubit system representing:

    S (Source): The black hole interior.

    H (Hidden): The hidden Hilbert domain.

    R (Radiation): The escaping Hawking radiation.

Unlike an evolution simulation, this script models a long-term evaporation process by calculating the entropy/information for a single emission step and then cumulating (summing) these identical, independent values over many steps (Nsteps).

The primary metrics it plots are:

    Cumulative Radiation Entropy (SR​(k)): The total entropy of all radiation collected up to step k.

    Cumulative Mutual Information (I(R:S∪H)): The total mutual information shared between the collected radiation and the combined Source/Hidden domains.

This script is responsible for generating the "Toy Page-like Curves" seen in the research.

Core Concepts Modeled

    build_step_density(mu): Creates the quantum state for a single emission. It models the core hypothesis: the information is in a superposition of being in the Source-Radiation subsystem (with probability 1−μ) and the Hidden-Radiation subsystem (with probability μ).

    dephasing_on_R(rho, p): Applies a noise channel to the Radiation qubit, modeling observation.

    step_entropies(mu, p): Calculates the von Neumann entropy and mutual information for a single emission step.

    simulate(Nsteps, mu, p): The main loop that calls step_entropies repeatedly and calculates the cumulative sum of the resulting values over Nsteps.

Requirements

The script requires the following Python libraries:

    numpy (for quantum state linear algebra)

    matplotlib (for generating the plots)

You can install the dependencies using pip:
pip install numpy matplotlib

How to Run

You can run the script directly from your terminal.
python3 quantum_entanglement_information_sim.py

Outputs

This script will open interactive matplotlib windows to display the plots. It does not save the plots to files automatically.

It generates four plots in total:

    Cumulative Radiation Entropy vs. Emission steps (k) for p=0.0

    Cumulative Radiation Entropy vs. Emission steps (k) for p=0.5

    Cumulative Mutual Information vs. Emission steps (k) for p=0.0

    Cumulative Mutual Information vs. Emission steps (k) for p=0.5

Each plot shows four lines corresponding to the different migration ratios (mu = 0.0, 0.3, 0.6, 0.9).
