


# [Input](@id man-input)

The file `src/Input.jl` and files contained in `src/input/` are collectively
referred to as the input module.
Here, all calculations prior to the finite element analysis are carried
out—including
- specification of the problem scenario, choice of mesh motion, and other
  parameters
- generation of the finite element mesh
- calculation of the basis functions and their derivatives
- initialization of the surface position and degrees of freedom
- application of boundary conditions

All files below are in the `src/input/` directory.


## [Enums.jl](@id man-enums)

The following enums are defined to simplify problem input.

```@docs
Scenario
```

```@docs
Topology
```

```@docs
Motion
```

```@docs
Boundary
```

```@docs
Corner
```

```@docs
Neumann
```

```@docs
Curve
```



## [Params.jl](@id man-params)

The requisite parameters for each simulation are split between the
[`Params`](@ref) struct and the keyword arguments `args`.
The data contained in [`Params`](@ref) is meant to be immutable, even across
restarted (or continued) simulations.
For this reason, temporal data (the initial time `t0`, initial time ID `t0_id`,
and time steps `Δts`) is treated as a keyword argument, despite it being
required for every simulation.
Additional details about required keywords are contained in
[`check_params`](@ref MembraneAleFem.check_params).

```@docs
MembraneAleFem.Params
```

```@docs
MembraneAleFem.check_params
```



## [Dof.jl](@id man-module-dof)

Within our ALE formalism, the membrane position ``\boldsymbol{x}`` is not a
fundamental unknown.
The mesh velocity ``\boldsymbol{v}^{\text{m}}``, on the other hand, *is* a
fundamental unknown, and over a time step ``\Delta t`` the membrane position is
updated according to
```math
\boldsymbol{x} (\zeta^\alpha, t + \Delta t)
\, = \, \boldsymbol{x} (\zeta^\alpha, t)
\, + \, \Delta t \, \boldsymbol{v}^{\text{m}} (\zeta^\alpha, t + \Delta t)
~,
```
where ``\Delta t`` is one of the entries of the array `ΔTS` specified in
[`Params.jl`](@ref man-params).
The membrane position is thus treated differently from the fundamental unknowns.
For this reason, we define all *possible* fundamental unknowns and the three
Cartesian components of the position in the `Dof` module:


```@docs
Dof.Unknown
```

```@docs
Dof.Position
```



## [Spline.jl](@id man-spline)

All B-spline calculations independent from finite element analysis, based
entirely on [**The NURBS Book**](https://doi.org/10.1007/978-3-642-97385-7)
[piegl-tiller](@citep).
See in particular *Chapter Two* for a description of B-spline functions, and
*Chapter Three* for an explanation of how the basis functions are used to
generate 1D curves and 2D surfaces in ``\mathbb{R}^3``.
We employ the same terminology in our code.

```@docs
KnotVector
```

```@docs
MembraneAleFem.knot_vector
```

```@docs
MembraneAleFem.get_knot_span_index
```

```@docs
MembraneAleFem.get_bspline_vals
```

```@docs
MembraneAleFem.get_bspline_ders
```

```@docs
MembraneAleFem.get_bspline_indices
```

```@docs
MembraneAleFem.get_1d_bspline_cps
```

```@docs
MembraneAleFem.get_2d_bspline_cps
```

```@docs
MembraneAleFem.collocate_ζ
```

```@docs
MembraneAleFem.get_unique_1d_elements
```


## [GaussPoint.jl](@id man-gauss-point)

As is standard in finite element analysis, Gauss--Legendre quadrature is used to
approximate integrals over the parametric domain.
Here, the primary objective is to calculate the 1-D integral
```math
I_1
\, = \, \int_{ζ_{\texttt{lo}}}^{ζ_{\texttt{hi}}} \!\!\! f(ζ) \, \text{d} ζ
~.
```
The domain of integration is mapped from ``[ ζ_{\texttt{lo}}, ζ_{\texttt{hi}} ]``
to ``[-1, +1]`` via the change of variables
```math
ζ
\, = \, \dfrac{1}{2} \, ξ \, \big(
	ζ_{\texttt{hi}}
	- ζ_{\texttt{lo}}
\big)
\, + \, \dfrac{1}{2} \, \big(
	ζ_{\texttt{hi}}
	+ ζ_{\texttt{lo}}
\big)
~,
```
for which the integral ``I_1`` is equivalently given by
```math
I_1
\, = \, \dfrac{1}{2} \, \big(
	ζ_{\texttt{hi}}
	- ζ_{\texttt{lo}}
\big) \int_{-1}^{+1} \!\! f \bigg(
	\dfrac{1}{2} \Big[
		ξ \big(
			ζ_{\texttt{hi}}
			- ζ_{\texttt{lo}}
		\big)
		+ \big(
			ζ_{\texttt{hi}}
			+ ζ_{\texttt{lo}}
		\big)
	\Big]
\bigg) \, \text{d} ξ
\, = \, \dfrac{1}{2} \, \big(
	ζ_{\texttt{hi}}
	- ζ_{\texttt{lo}}
\big) \int_{-1}^{+1} \!\! \hat{f} (ξ) \, \text{d} ξ
~.
```
At this point, the integral is approximated by applying the Gaussian quadrature
formula
```math
I_1
\, ≈ \, \dfrac{1}{2} \, \big(
	ζ_{\texttt{hi}}
	- ζ_{\texttt{lo}}
\big) \, \sum_{k = 1}^{\texttt{ngp}} \hat{w}_k \, \hat{f} (ξ_k)
~,
```
where ``\texttt{ngp}`` is the number of 1D Gauss points, ``\hat{w}_k`` are the
corresponding weights, and ``ξ_k`` are the associated Gauss--Legendre points.
In our code, [`GaussPointsξ`](@ref) contains ``\texttt{ngp}``, ``\hat{w}_k``,
and ``ξ_k``.
Note that this `struct` is the same for every line element, because 
``ξ ∈ [-1, +1]``
always.


```@docs
GaussPointsξ
```


In our finite element implementation, it is often convenient to evaluate
integrals over the original variable `ζ` rather than the transformed variable
`ξ`.
To this end, we recognize that
```math
I_1
\, ≈ \, \sum_{k = 1}^{\texttt{ngp}} w_k \, f (ζ_k)
~,
```
where
```math
w_k
\, := \,\dfrac{1}{2} \, \hat{w}_k \, \big(
	ζ_{\texttt{hi}}
	- ζ_{\texttt{lo}}
\big)
\qquad
\text{and}
\qquad
ζ_k
\, := \, \dfrac{1}{2} \Big[
	ξ_k \big(
		ζ_{\texttt{hi}}
		- ζ_{\texttt{lo}}
	\big)
	+ \big(
		ζ_{\texttt{hi}}
		+ ζ_{\texttt{lo}}
	\big)
\Big]
~.
```
In our code, [`GaussPointsζ`](@ref) contains `ngp`, ``w_k``, and ``ζ_k``.

```@docs
GaussPointsζ
```


## [GpBasisFn.jl](@id man-gp-basis-fn)

Calculation of nonzero B-spline basis functions at every Gaussian quadrature
point.
Since the basis functions are determined upon discretizing the parametric
domain, they are calculated once and then stored in the structs listed
below—which systematically build up in complexity.

```@docs
GpBasisFnsζ
```

```@docs
GpBasisFnsζα
```

```@docs
LineGpBasisFns
```

```@docs
BdryGpBasisFns
```

```@docs
AreaGpBasisFns
```


## [Mesh.jl and Bc.jl](@id man-mesh)

The mesh is one of the main constructions required for finite element analysis.
In this codebase, the functions associated with mesh generation are divided into
two files.
`Mesh.jl` deals with aspects of mesh generation that are independent of the
scenario under consideration, while `Bc.jl` handles the boundary conditions and
their effect on mesh organization—which is scenario-dependent.


### Mesh.jl

All relevant information is contained in the [`Mesh`](@ref) struct, and will not
change over the course of a simulation.

```@docs
Mesh
```

```@docs
MembraneAleFem.generate_scenario
```

The following helper functions allow the user to easily obtain basis functions,
their derivatives, and Gauss point weights both in the mesh interior and on the
mesh boundary.

```@docs
MembraneAleFem.get_basis_fns
```

```@docs
MembraneAleFem.get_gpw
```

```@docs
MembraneAleFem.get_N
```

```@docs
MembraneAleFem.get_∂Nα
```

```@docs
MembraneAleFem.get_∂∂Nαβ
```



### [Bc.jl](@id man-bc)

When constructing the [`Mesh`](@ref),
[`generate_scenario`](@ref MembraneAleFem.generate_scenario) (in `Mesh.jl`)
calls [`get_scenario_bc_info`](@ref)—which is a helper function that then calls
the appropriate [`Scenario`](@ref)-specific function.
Each such function returns the following:
- `dofs` --> global ordering of active [`Unknown`](@ref Dof.Unknown)s
- `ndf`  --> number of active [`Unknown`](@ref Dof.Unknown)s
- `ID`   --> matrix mapping node number to indices of unknowns at that node
- `inh_dir_bcs` --> list of inhomogeneous Dirichlet boundary conditions
- `inh_neu_bcs` --> list of inhomogeneous Neumann boundary conditions


```@docs
get_scenario_bc_info
```

```@docs
MembraneAleFem.get_f_cavi_bc_info
```

```@docs
MembraneAleFem.get_f_coue_bc_info
```

```@docs
MembraneAleFem.get_f_pois_bc_info
```

```@docs
MembraneAleFem.get_f_pull_bc_info
```

```@docs
MembraneAleFem.get_f_bend_bc_info
```

```@docs
MembraneAleFem.get_dofs
```



## [Input.jl](@id man-input-jl)


```@docs
MembraneAleFem.prepare_input
```

