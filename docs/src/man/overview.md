# [Overview](@id man-overview)

In this manual, all paths are provided relative to the root directory
`MembraneAleFem.jl/`.
The source code is provided in `src/`, which contains the following files and
folders:
```code
$ tree src/
src/
├── analysis
│   ├── FiniteElement.jl
│   ├── GeoDynStress.jl
│   └── PullForce.jl
├── Analysis.jl
├── input
│   ├── Bc.jl
│   ├── Dof.jl
│   ├── Enums.jl
│   ├── GaussPoint.jl
│   ├── GpBasisFn.jl
│   ├── Mesh.jl
│   ├── Params.jl
│   └── Spline.jl
├── Input.jl
├── MembraneAleFem.jl
└── Output.jl
```

The function [`solve`](@ref) in the main program file `src/MembraneAleFem.jl`
carries out all calculations, which are divided into three parts:
1. [Input](@ref man-input) (pre-processing)
2. [Analysis](@ref man-analysis) (solution)
3. [Output](@ref man-output) (post-processing)
We briefly overview the three modules before presenting further details below.

`src/Input.jl` and corresponding files in `src/input/` deal with the problem
set-up, mesh generation, and degree-of-freedom initialization.
Since the parametrization ``\zeta^\alpha`` of the parametric domain
``\Omega \subset \mathbb{R}^2``
does not change, all basis functions are determined once ``\Omega`` is initially
discretized.
We employ B-spline basis functions, as they are well-characterized and satisfy
``C^1``-continuity across finite elements.
The local-to-global mapping between elements and degrees of freedom is also
determined here.

All finite element analysis is carried out in `src/Analysis.jl` and the files in
the `src/analysis/` folder.
In particular, `src/analysis/GeoDynStress.jl` takes as arguments the elemental
degrees of freedom, and calculates all relevant geometric and dynamic
quantities—including the stresses and couple-stresses of the membrane and mesh.
The file `src/analysis/FiniteElement.jl` carries out all finite element
analysis, including the calculation of the residual vector and tangent diffusion
matrix.
`src/analysis/PullForce.jl` calculates the pull force on a tether, and is thus
specific to one of the available scenarios.

Currently, `src/Output.jl` includes a minimal set of functions that process the
results of our simulations.



## [MembraneAleFem.jl](@id man-membrane-ale-fem)

The main file `src/MembraneAleFem.jl` contains a single function, `solve()`,
which then calls relevant functions in the aforementioned files.

```@docs
solve
```

