module Dof


"""
    Dof.Unknown

Possible degrees of freedom (or fundamental unknowns) at each node.

For different mesh motions, different fundamental unknowns are required.
A Lagrangian simulation requires only `vx`, `vy`, `vz`, and `λ`, while an ALE
simulation requires all eight fundamental unknowns.
Currently, the mapping from the choice of [`Motion`](@ref MembraneAleFem.Motion)
to the corresponding fundamental unknowns is specified in the function
[`get_dofs`](@ref MembraneAleFem.get_dofs) and stored in
[`Mesh`](@ref MembraneAleFem.Mesh)`.dofs`.

This enum is encapsulated in a module; access elements with e.g. `Dof.vx`.

Types of fundamental unknowns:
- `vx`  --> ``x``-velocity,      ``v_x``
- `vy`  --> ``y``-velocity,      ``v_y``
- `vz`  --> ``z``-velocity,      ``v_z``
- `vmx` --> ``x``-mesh velocity, ``v^{\\text{m}}_x``
- `vmy` --> ``y``-mesh velocity, ``v^{\\text{m}}_y``
- `vmz` --> ``z``-mesh velocity, ``v^{\\text{m}}_z``
- `λ`   --> surface tension,     ``λ``
- `pm`  --> mesh pressure,       ``p^{\\text{m}}``
"""
@enum Unknown begin
  vx  = 1; vy  = 2; vz  = 3;
  vmx = 4; vmy = 5; vmz = 6;
  λ   = 7;
  pm  = 8;
end # Unknown


"""
    Dof.Position

Cartesian components of the mesh position at each node.

This enum is encapsulated in a module; access elements with e.g. `Dof.xm`.

Options:
- `xm` --> ``x``-position
- `ym` --> ``y``-position
- `zm` --> ``z``-position
"""
@enum Position begin
  xm  = 1; ym  = 2; zm  = 3;
end # Position

end # module Dof
