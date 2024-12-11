
# MembraneAleFem.jl

Source code for arbitrary Lagrangian--Eulerian (ALE) lipid membrane finite
element method (FEM).


## Overview

**MembraneAleFem.jl** solves the continuum equations governing the dynamics of a
lipid membrane, using an arbitrary Lagrangian--Eulerian (ALE) finite element
method (FEM).
More details, including the entire codebase, will soon be pushed to this
repository.
In the meantime, please see the description of the method in the
[arXiv preprint](https://arxiv.org/pdf/2412.07596).


## Results

One of the capabilities of the codebase is to pull a tether from an initially
flat membrane patch.
The following three videos capture the results of Lagrangian, Eulerian, and ALE
simulations---corresponding to Figures 4 and 5 of the manuscript.


### Lagrangian tether pulling

https://github.com/user-attachments/assets/5b345ea8-f9b6-44b6-b673-c7658df16aa8


### Eulerian tether pulling

https://github.com/user-attachments/assets/578ed295-59ab-43ed-9e98-a2f57a590250


### ALE tether pulling

https://github.com/user-attachments/assets/abdb39f2-d9be-4f8d-acec-97cff811b1b8


### ALE tether translation

For the first time, one can translate the tether laterally across the membrane
surface.
Doing so is only possible with an ALE mesh motion, and is shown in the following
video.

https://github.com/user-attachments/assets/ed068922-7301-4c63-8ec0-5ba4ee935e0d

