

"""
    Scenario

Scenarios to simulate: `instances(Scenario)`

- `F_CAVI = 1` --> lid-driven cavity
- `F_COUE = 2` --> Couette flow
- `F_POIS = 3` --> Hagen--Poiseuille flow
- `F_PULL = 4` --> tether pulling from a flat patch
- `F_BEND = 5` --> bending of a flat patch

The first letter indicates mesh [`Topology`](@ref): `F` for a flat patch and `C` for a
cylinder.
At present, no cylindrical scenarios are implemented.

To create a new scenario:
- add and export a new entry in the `Scenario` enum
- edit [`get_scenario_bc_info`](@ref) in `src/input/Bc.jl`
- write the corresponding `apply_<SCENARIO>_bcs` function in `src/input/Mesh.jl` 
- add the new `scenario` to the `implemented` list in `src/input/Params.jl`
"""
@enum Scenario begin
  F_CAVI = 1;
  F_COUE = 2;
  F_POIS = 3;
  F_PULL = 4;
  F_BEND = 5;
end

export Scenario, F_CAVI, F_COUE, F_POIS, F_PULL, F_BEND


"""
    Topology

Mesh topologies: `instances(Topology)`

- `FLAT     = 1`
- `CYLINDER = 2`
"""
@enum Topology begin
  FLAT     = 1;
  CYLINDER = 2;
end

export Topology, FLAT, CYLINDER


"""
    Motion

Types of mesh motion: `instances(Motion)`

- `STATIC = 1` --> mesh does not move
- `EUL    = 2` --> Eulerian mesh motion
- `LAG    = 3` --> Lagrangian mesh motion
- `ALEV   = 4` --> ALE-viscous mesh motion
- `ALEVB  = 5` --> ALE-viscous-bending mesh motion

The choice of motion is stored in `Params.motion` variable in
[`Params.jl`](@ref man-params), and is relevant both to the specification of
boundary conditions in [`Bc.jl`](@ref man-mesh) and the
[finite element analysis](@ref man-analysis).
"""
@enum Motion begin
  STATIC = 1;
  EUL    = 2;
  LAG    = 3;
  ALEV   = 4;
  ALEVB  = 5;
end

export Motion, STATIC, EUL, LAG, ALEV, ALEVB


"""
    Boundary

Enumerate four mesh boundaries: `instances(Boundary)`, never to be used as an
index.

- `BOTTOM = 1`
- `RIGHT  = 2`
- `TOP    = 3`
- `LEFT   = 4`
"""
@enum Boundary begin
  BOTTOM = 1;
  RIGHT  = 2;
  TOP    = 3;
  LEFT   = 4;
end

export Boundary, BOTTOM, RIGHT, TOP, LEFT


"""
    Corner

Enumerate four mesh corners: `instances(Corner)`, never to be used as an index.

- `BOTTOM_LEFT  = 1`
- `BOTTOM_RIGHT = 2`
- `TOP_LEFT     = 3`
- `TOP_RIGHT    = 4`
"""
@enum Corner begin
  BOTTOM_LEFT  = 1;
  BOTTOM_RIGHT = 2;
  TOP_LEFT     = 3;
  TOP_RIGHT    = 4;
end

export Corner, BOTTOM_LEFT, BOTTOM_RIGHT, TOP_LEFT, TOP_RIGHT


"""
    Neumann

Enumerate types of Neumann boundary conditions: `instances(Neumann)`

- `SHEAR     = 1` --> apply in-plane force/length in the τ direction
- `STRETCH   = 2` --> apply in-plane force/length in the ν direction
- `MOMENT    = 3` --> apply boundary moment

Associate finite element calculations are found in
[`calc_bdry_element_residual`](@ref)
"""
@enum Neumann begin
  SHEAR     = 1;
  STRETCH   = 2;
  MOMENT    = 3;
end

export Neumann, SHEAR, STRETCH, MOMENT


"""
    Curve

Enumerate types of curves: `instances(Curve)`

- `CLAMPED`
- `CLOSED`

These curves are used to generate B-spline basis functions—see
[`Spline.jl`](@ref man-spline)
"""
@enum Curve begin
  CLAMPED = 1;
  CLOSED  = 2;
end

export Curve, CLAMPED, CLOSED

