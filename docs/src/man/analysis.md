


# [Analysis](@id man-analysis)

The file `src/Analysis.jl` and files contained in the `src/analysis/` directory
are referred to as the analysis module.
All finite element analysis is carried out here, including
- determination of the membrane geometry and dynamics
- calculation of the residual vector and tangent matrix
- implementation of Newton--Raphson iteration, and time stepping
- calculation of the pull force


## Overview

Though several different mesh motions are implemented, we provide details for
only the most involved one: the ``\texttt{ALE-vb}`` mesh motion.
At any time step, we begin with the vector of known degrees of freedom
``[\mathbf{u} (t)]``, and seek to calculate the unknown degrees of freedom
``[\mathbf{u} (t + \Delta t)]``.
In practice, we work primarily with the degrees of freedom corresponding to the
various unknowns—in this case, the velocities ``[\mathbf{v} (t)]``, mesh
velocities ``[\mathbf{v}^{\text{m}} (t)]``, surface tensions
``[\bm{\lambda} (t)]``, and mesh pressures ``[\mathbf{p}^{\text{m}} (t)]``.

The direct Galerkin expression for the ``\texttt{ALE-vb}`` mesh motion is
provided in §2.7.3 of our manuscript [sahu-arxiv-2024](@citep).
Upon discretization, the direct Galerkin expression is satisfied by solving the
set of equations
```math
[\mathbf{r}^\lambda]
\, = \, [\mathbf{0}]
~,
\qquad
[\mathbf{r}^{\text{v}}]
\, = \, [\mathbf{0}]
~,
\qquad
[\mathbf{r}^{\text{m}}]
\, = \, [\mathbf{0}]
~,
\qquad
\text{and}
\qquad
[\mathbf{r}^{\text{p}}]
\, = \, [\mathbf{0}]
~,
```
where each vector corresponds to the portion of the residual vector
``[\mathbf{r}]`` associated with a particular unknown.
The unknown degree of freedom vector ``[\mathbf{u} (t + \Delta t)]`` then
satisfies
```math
[\mathbf{r} ( [\mathbf{u} (t + \Delta t)] ) ]
\, = \, [\mathbf{0}]
~.
```
The above equation is solved with the Newton--Raphson method, in which our
initial approximation for ``[\mathbf{u} (t + \Delta t)]``, denoted
``[\mathbf{u} (t + \Delta t)]_0``, is chosen to be ``[\mathbf{u} (t)]``.
Successive approximations are determined according to
```math
[\mathbf{u} (t + \Delta t)]_{j+1}
\, = \, [\mathbf{u} (t + \Delta t)]_j
\, - \, [\mathbf{K}]_j^{-1} [\mathbf{r} ( [\mathbf{u} (t + \Delta t)]_j ) ]
~,
```
where the global tangent diffusion matrix at the ``j^{\text{th}}`` iteration is
given by
```math
[\mathbf{K}]_j
\, := \, \dfrac{\partial [\mathbf{r}]}{\partial [\mathbf{u}]}
\bigg\rvert_{[\mathbf{u} (t + \Delta t)]_j}
~.
```
Note that the diffusion matrix is calculated numerically by extending
``[\mathbf{r}]`` and ``[\mathbf{u}]`` into the complex plane
\[[lyness-siam-1967](@citet), [lyness-mc-1968](@citet)\].
Details of each portion of our calculation can be found below.


## [GeoDynStress.jl](@id man-geo-dyn-stress)

To calculate the residual vector, integrals over the parametric domain are
evaluated as the sum of integrals over each finite element.
The latter are numerically calculated via Gauss point integration, for which
(see [GaussPoint.jl](@ref man-gauss-point) and
[GpBasisFn.jl](@ref man-gp-basis-fn))
```math
\int_\Omega f(\zeta^\alpha) ~ \text{d}\Omega
\ = \ \sum_{e=1}^{\texttt{nel}} \int_{\Omega^e} f(\zeta^\alpha) ~ \text{d}\Omega
\ = \ \sum_{e=1}^{\texttt{nel}} \sum_{k=1}^{\texttt{ngp}}
	w_k f(\zeta^{\alpha, e}_k)
~.
```
To simplify such calculations, the struct [`GeoDynStress`](@ref) contains all
quantities at a single Gauss point that are used in the determination of the
residual vector.
Many of the quantities are self-explanatory; here we describe the few that are
not.

To begin, we use Voigt notation to represent `2×2` matrices as `3×1` vectors,
where we take advantage of symmetries in the stresses, couple stresses, and
arbitrary variations.
Thus,

- `𝞂` is a `3×1` vector containing ``(σ^{1 1}, σ^{2 2}, 2 σ^{1 2})``


- `𝞂m` is a `3×1` vector containing
  ``(σ^{1 1}_{\text{m}}, σ^{2 2}_{\text{m}}, 2 σ^{1 2}_{\text{m}})``


- `𝗠` is a `3×1` vector containing ``(M^{1 1}, M^{2 2}, M^{1 2} + M^{2 1})``


In addition, we introduce the matrices `𝗕a` and `𝗕b` containing products of the
basis functions with surface geometry.
More precisely, we have
```math
[\mathbf{B}^a]
\, = \, \bigg[ ~
	[\mathbf{B}^a_1] ~
	[\mathbf{B}^a_2] ~
	\ldots ~
	[\mathbf{B}^a_{\text{nen}}] ~
\bigg]
\qquad
\text{and}
\qquad
[\mathbf{B}^b]
\, = \, \bigg[ ~
	[\mathbf{B}^b_1] ~
	[\mathbf{B}^b_2] ~
	\ldots ~
	[\mathbf{B}^b_{\text{nen}}] ~
\bigg]
~,
```
where
```math
[\mathbf{B}^a_j]
\, := \, \begin{bmatrix}
	[\bm{a}_1]^{\text{T}} N_{j, 1} \\[4pt]
	[\bm{a}_2]^{\text{T}} N_{j, 2} \\[4pt]
	(
		[\bm{a}_1]^{\text{T}} N_{j, 2}
		\, + \, [\bm{a}_2]^{\text{T}} N_{j, 1}
	) / 2 \\[2pt]
\end{bmatrix}
\qquad
\text{and}
\qquad
[\mathbf{B}^b_j]
\, := \, \begin{bmatrix}
	[\bm{n}]^{\text{T}} N_{j; 1 1} \\[4pt]
	[\bm{n}]^{\text{T}} N_{j; 2 2} \\[4pt]
	[\bm{n}]^{\text{T}} \, (
		N_{j; 1 2}
		\, + \, N_{j; 2 1}
	) / 2 \\[2pt]
\end{bmatrix}
```
As we will see, determining these quantities greatly simplifies our later finite
element calculations.


```@docs
GeoDynStress
```

```@docs
MembraneAleFem.get_dof_cps
```


## [FiniteElement.jl](@id man-finite-element)

```@docs
MembraneAleFem.time_step!
```

```@docs
MembraneAleFem.calc_r_K
```

```@docs
MembraneAleFem.calc_elem_residual
```

```@docs
MembraneAleFem.calc_elem_dof_residuals
```

```@docs
MembraneAleFem.calc_bdry_element_residual
```

```@docs
MembraneAleFem.update_xms!
```

```@docs
MembraneAleFem.calc_τ_ν
```


## [PullForce.jl](@id man-pull-force)

Specific calculations related to determining the force during tether pulling, as
detailed in Appendix C of [sahu-arxiv-2024](@citet).

```@docs
MembraneAleFem.get_pull_el_id
```

```@docs
MembraneAleFem.get_adj_maps
```

```@docs
MembraneAleFem.calc_pull_force
```


## [Analysis.jl](@id man-analysis-jl)

```@docs
MembraneAleFem.run_analysis!
```



