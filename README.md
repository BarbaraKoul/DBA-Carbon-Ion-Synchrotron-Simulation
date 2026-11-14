# Simulation of Double-Bend Achromat Carbon-Ion Synchrotron for cancer therapy
<img width="640" height="480" alt="Figure_1" src="https://github.com/user-attachments/assets/33675a59-0543-4be0-af0a-66fecb335560" />

## Summary
Based on the paper of Xuanhao Zhang et al. [Lattice Design of a Carbon-Ion Synchrotron based on Double-Bend Achromat Lens](https://doi.org/10.48550/arXiv.2007.11787),
in this project we reproduce the proposed design working with [Xsuite](https://doi.org/10.18429/JACoW-HB2023-TUA2I1),  a collection of python packages for the simulation of the beam dynamics in particle accelerators developed in CERN. More specifically, we focus on:

- Plotting the ring with its elements (Dipoles, Quadrupoles, Sextupoles etc.).
- Calculating and plotting the Twiss parameters, the dispersion and the chromaticity of the beam.
- Searching the drift lengths between the elements with matching in order to reach the values of the paper.
- Matching the quadrupole strengths (kq<sub>f1</sub>, kq<sub>d1</sub>, kq<sub>f2</sub>, kq<sub>d2</sub>) and the sextupole for further optimisation.
- Optimisation of the synchrotron by benchmarking Simulated Annealing Algorithm (SAA), Genetic Algorithms (GA), such as NSGA-II and NSGA-III and Bayesian optimisation.
- Defining the RF Cavity traits (Voltage, Harmonic, Size, Placement etc.).
- Calculating and optimising the tune Q<sub>s</sub>.

## How to use
>[!IMPORTANT]
> Linux OS or WSL for Windows users is mandatory to run this project.

Before any installment, it is preferable to create a Virtual Environment for the project. We propose **Miniforge**(Conda-forge ecosystem).

Install **Xsuite** as stated in the [User's guide](https://xsuite.readthedocs.io/en/latest/usersguide.html) and make sure you have **Python3** version or later. Then, clone this repository:
```
bash git clone <https://github.com/BarbaraKoul/DBA-Carbon-Ion-Synchrotron-Simulation>
```
You can run the python files in your code editor or terminal. 
>[!NOTE]
> To see the diagrams, make sure you have **PyQt5** version or later. Otherwise, you can save the plots as **.png** files in your folder.

## Lattice Design
For the ring, we made 2 files: ``lattice_design.py`` and ``lattice_design_RF.py``. 

The ``lattice_design.py`` is created for:
- The creation of an environment and the declaration of the variables of the lattice (Magnets, Drift lengths, particle properties).
- The ring plot.
- The calculation and plotting of the Twiss parameters with the 4D method, the dispersion and the chromaticity.
- The calculation of the tunes Q<sub>x</sub>, Q<sub>y</sub>.

The ``lattice_design_RF.py`` is created for:

## Matching
For matching, we have again 2 files: ``matching.py`` and ``matching_strengths.py``.

The ``matching.py`` is created for estimating the distance between the elements, seeing that they are not stated in the publication:
- We exploit the 180<sup>o</sup> symmetry of the ring and we reduce the variable lengths to 7.
- Our targets in this case are: Q<sub>x</sub>, Q<sub>y</sub>, circumference, &beta;<sub>x</sub>, &beta;<sub>y</sub>, D<sub>x</sub> and &gamma;<sub>t</sub>.

The ``matching_strengths.py`` is created for optimisation of the Twiss parameters of the beam by varying Quadrupole strengths:
- Our targets in this case are: Q<sub>x</sub>, Q<sub>y</sub>, &beta;<sub>x</sub> and D<sub>x</sub>=0 at the straight sections of the accelerator.

## Simulated Annealing Algorithm for Quadrupole Strength Optimisation


# RF Cavities
