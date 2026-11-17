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
-The addition of an RF Cavity to the previously defined environment of the simple lattice.
-The calculation of the tune longitudinal to the particle, Q<sub>s</sub>.
-The plotting over time of RMS Beam Size on the Z axis and the Twiss parameters.

## Matching
For matching, we have again 2 files: ``matching.py`` and ``matching_strengths.py``.

The ``matching.py`` is created for estimating the distance between the elements, seeing that they are not stated in the publication:
- We exploit the 180<sup>o</sup> symmetry of the ring and we reduce the variable lengths to 7.
- Our targets in this case are: Q<sub>x</sub>, Q<sub>y</sub>, circumference, &beta;<sub>x</sub>, &beta;<sub>y</sub>, D<sub>x</sub> and &gamma;<sub>t</sub>.

The ``matching_strengths.py`` is created for optimisation of the Twiss parameters of the beam by varying Quadrupole strengths:
- Our targets in this case are: Q<sub>x</sub>, Q<sub>y</sub>, &beta;<sub>x</sub> and D<sub>x</sub>=0 at the straight sections of the accelerator.

## Simulated Annealing Algorithm for Quadrupole Strength Optimisation
As stated in the [Cornell University Wiki:](https://optimization.cbe.cornell.edu/index.php?title=Simulated_annealing)</br>
*SA is a probabilistic technique for finding approximate solutions to optimization problems, particularly those with large search spaces. Inspired by the annealing process in metallurgy, the algorithm explores the solution space by occasionally accepting worse solutions with a probability that diminishes over time. This reduces the risk of getting stuck in local optima and increases the likelihood of discovering the global optimum.*</br>

A pseudocode for this algorithm is the following:</br>

```
1: T = T0
2: Ts = T1
3: δ = δ0
4: x = x0
5: while T > Ts do
6:  x1 = RandSelect(neighbourhood(x))
7:  ∆F = Fitness(x) − Fitness(x1)
8:  if ∆F > 0 then
9:    x ← x1
10: else if random(0,1) < exp(−∆F/T) then
11:   x ← x1
12: T = T − δ
```
## RF Cavities
In order to avoid interfering with the 4D Twiss methods, two new files were created for the RF Cavity simulation (or 6D Twiss method). ``lattice_design_RF.py`` and ``configRF.py``. </br>
In ``configRF.py``: </br>
* We define the the RF Voltage as 55kV </br>
* The RF frequency from a calculation that was made on hand, possibly on a wrong assumption. </br> 
* and the RF harmonic as 32, which is a way to calculate the required frequency through Xsuite </br>
We also defined the RF Cavity length and the rest is left as is from ``config.py`` </br>

In ``lattice_design_RF.py`` the largest change being made is of course the Twiss 6D method. The ring is mostly similar, except for the fact that we slice it to 
insert 'rf1' or the RF cavity in the lattice 1.874 meters left from the sextupole. The tune is no longer 0 so we print it in this mode. 6D method allows us to see
how many of the values we considered stable before, change with time. BeamStatMonitor allows us to do that. We define an infinitesimal element of the lattice that
allows us to take measurements every time a particle (or rather some bunched particles) pass through that element, which we will place at 0. </br>

We first however have to define the injection energy, normalized emittance on x and y as well as rms beam size on z axis and finally the number of particles. 
Xsuite runs out of memory fast with the number of particles however. 10^5 particles would need 74.5GB of RAM to run, we therefore keep our particle numbers low,
around the low hundreds to not run into any problems. The injection energy can be defined by the max penetration depth together with the realistic capabilities
of an RF cavity on a 55 meter synchrotron. Normalized emittance both on x and y and rms beam size on z can be defined AFTER evaluating the results, in other words
it is a better approach to use the results to define the starting values than the other way around. After all emittance on injection is defined by the injection
linac which we do not have data on (same with rms beam size). </br>

We then plot over every turn:</br>
the rms beam size on z axis, the deviation of the beam from zeta, the betatron x function, the betatron y function, the dispersion on X axis and the dispersion on 
Y axis respectively. </br> 
