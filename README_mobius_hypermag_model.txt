
Möbius Hypermagnetic Model

mobius_hypermag_model.py

Dr. Chris Caulkins
c.caulkins@FreeRangeThinker.org

Overview

This Python script implements a numerical simulation of the Möbius-Hilbert framework, a theoretical model for resolving the black hole information paradox. It simulates a 3-qubit quantum system representing:

    S (Source): The black hole interior.

    H (Hidden): The hidden Hilbert domain.

    R (Radiation): The escaping Hawking radiation.

The script's primary purpose is to model the flow of quantum information (measured by mutual information and von Neumann entropy) between these three subsystems under the influence of:

    A "hypermagnetic" rotation acting on the radiation.

    A "Möbius transfer" (an iSWAP-like entangling gate) between the hidden domain and the radiation.

    A dephasing (observation) channel acting on the radiation.

The script is non-interactive and is designed to run, generate three specific plots that demonstrate the model's key principles, and then save them as PNG files in the same directory.

Requirements

The script requires the following Python libraries:

    numpy (for quantum state linear algebra)

    matplotlib (for generating the plots)

You can install the dependencies using pip:
Bash

pip install numpy matplotlib

How to Run

You can run the script directly from your terminal.
Bash

python3 mobius_hypermag_model.py

Note: This script uses the matplotlib.use("Agg") backend. This means it will not open any GUI windows to display the plots. Instead, it will render the plots directly to files.

Outputs

When executed, the script will run three simulations and save the following three files in the same directory:

    MI_RH_RS_vs_muB.png:

        Description: A plot of the mutual information I(R:H) vs. I(R:S) as the hypermagnetic rotation angle (μB​) is varied.

        Purpose: This demonstrates how the hypermagnetic field organizes information, showing that the radiation (R) becomes more correlated with the hidden domain (H) than the source (S) under certain conditions.

    MI_RH_RS_vs_p_for_muB.png:

        Description: A plot of I(R:H) vs. I(R:S) as the dephasing (observation) strength (p) is varied, at a fixed rotation angle.

        Purpose: This demonstrates the model's robustness to observational decoherence.

    entropy_R_page_like.png:

        Description: A plot of the von Neumann entropy of the radiation subsystem, S(R), as it evolves over a number of "toy" evaporation steps.

        Purpose: This demonstrates that the model's dynamics can reproduce a "Page-like curve," where the entropy of the radiation behaves in a way consistent with information conservation.

After saving the files, the script will print the full, absolute paths to the three generated images for easy access.