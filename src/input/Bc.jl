

export get_scenario_bc_info


"""
    get_scenario_bc_info(numnp, IX, bdry_nodes, ..., p::Params; args...)

Return [`Scenario`](@ref)-specific information for the provided `Params`.

Each of the functions called returns
`(dofs, ndf, ID, inh_dir_bcs, inh_neu_bcs)`,
in that order.
"""
function get_scenario_bc_info(
    numnp::Int64,
    IX::Matrix{Int64},
    bdry_nodes::Dict{Boundary, Vector{Int64}},
    bdry_inner_nodes::Dict{Boundary, Vector{Int64}},
    crnr_nodes::Dict{Corner, Int64},
    p::Params;
    args...
  )

  if p.scenario == F_CAVI
    return get_f_cavi_bc_info(numnp, bdry_nodes, crnr_nodes);
  elseif p.scenario == F_COUE
    return get_f_coue_bc_info(numnp, bdry_nodes);
  elseif p.scenario == F_POIS
    return get_f_pois_bc_info(numnp, bdry_nodes);
  elseif p.scenario == F_BEND
    return get_f_bend_bc_info(numnp, bdry_nodes, p; args...);
  elseif p.scenario == F_PULL
    return get_f_pull_bc_info(numnp, IX, bdry_nodes,
                              bdry_inner_nodes, p; args...);
  else
    @assert false "Need boundary conditions for $(scenario) scenario";
    return Nothing;
  end
end # get_scenario_bc_info


"""
    get_f_cavi_bc_info(numnp, bdry_nodes, crnr_nodes)

Generate [`Mesh`](@ref) data for the flat cavity [`Scenario`](@ref) `F_CAVI`.

In this scenario, a `STATIC` mesh [`Motion`](@ref) is required; relevant
[`Unknown`](@ref Dof.Unknown)s are `vx`, `vy`, and `λ`.
The following boundary conditions are prescribed:
- ``v_x = 1.0`` along the top edge (corners excluded)
- ``v_x = 0.0`` along all other edges (corners included)
- ``v_y = 0.0`` along all edges
- ``λ = 0.0`` at (or near) the center of the domain
Note that ``λ`` is pinned to remove the global indeterminacy of the surface
tension, up to a constant.
"""
function get_f_cavi_bc_info(
    numnp::Int64,
    bdry_nodes::Dict{Boundary, Vector{Int64}},
    crnr_nodes::Dict{Corner, Int64}
  )

  # fundamental unknowns
  dofs = Dict(Dof.vx=>1, Dof.vy=>2, Dof.λ=>3);
  ndf = length(values(dofs));

  ## construct ID matrix
  ID = zeros(Int64, ndf, numnp);

  ## remove Dirichlet boundary conditions
  inh_dir_bcs = Tuple{Dof.Unknown, Int64, Float64}[];

  for i in bdry_nodes[BOTTOM] # bottom edge
    ID[dofs[Dof.vx], i] = -1; # vx = 0
    ID[dofs[Dof.vy], i] = -1; # vy = 0
  end
  for i in bdry_nodes[TOP]    # top edge
    ID[dofs[Dof.vx], i] = -1; # vx = 1
    ID[dofs[Dof.vy], i] = -1; # vy = 0
    # exclude top corners for non-homogeneous velocity boundary condition
    if i != crnr_nodes[TOP_LEFT] && i != crnr_nodes[TOP_RIGHT]
      push!(inh_dir_bcs, (Dof.vx, i, 1.0)); # vx = 1
    end
  end
  # left edge
  for i in bdry_nodes[LEFT]   # left edge
    ID[dofs[Dof.vx], i] = -1; # vx = 0
    ID[dofs[Dof.vy], i] = -1; # vy = 0
  end
  # right edge
  for i in bdry_nodes[RIGHT]  # right edge
    ID[dofs[Dof.vx], i] = -1; # vx = 0
    ID[dofs[Dof.vy], i] = -1; # vy = 0
  end

  # pin surface tension in the center of the domain
  center_node_id = floor(Int64, numnp/2)+1;
  ID[dofs[Dof.λ], center_node_id] = -1;
  push!(inh_dir_bcs, (Dof.λ, center_node_id, 0.0));

  # Neumann boundary conditions
  inh_neu_bcs = Tuple{Boundary, Neumann, Float64}[];

  return dofs, ndf, ID, inh_dir_bcs, inh_neu_bcs;
end # get_f_cavi_bc_info


"""
    get_f_coue_bc_info(numnp, bdry_nodes)

Generate [`Mesh`](@ref) data for the flat Couette [`Scenario`](@ref) `F_COUE`.

In this scenario, a `STATIC` mesh [`Motion`](@ref) is required; relevant
[`Unknown`](@ref Dof.Unknown)s are `vx`, `vy`, and `λ`.
The following boundary conditions are prescribed:
- ``v_x = 0.0`` along the bottom edge (corners included)
- ``v_x = 3.0`` along the top edge (corners included)
- ``v_y = 0.0`` along all edges
- ``\\bm{\\bar{F}} = 4\\bm{\\nu}`` on the left and right edges
"""
function get_f_coue_bc_info(
    numnp::Int64,
    bdry_nodes::Dict{Boundary, Vector{Int64}}
  )

  # fundamental unknowns
  dofs = Dict(Dof.vx=>1, Dof.vy=>2, Dof.λ=>3);
  ndf = length(values(dofs));

  ## construct ID matrix
  ID = zeros(Int64, ndf, numnp);

  ## remove Dirichlet boundary conditions
  inh_dir_bcs = Tuple{Dof.Unknown, Int64, Float64}[];

  for i in bdry_nodes[BOTTOM] # bottom edge
    ID[dofs[Dof.vx], i] = -1; # vx = 0
    ID[dofs[Dof.vy], i] = -1; # vy = 0
  end
  for i in bdry_nodes[TOP]    # top edge
    ID[dofs[Dof.vx], i] = -1; # vx = 3
    ID[dofs[Dof.vy], i] = -1; # vy = 0
    push!(inh_dir_bcs, (Dof.vx, i, 3.0)); # vx = 3
  end
  for i in bdry_nodes[LEFT]   # left edge
    ID[dofs[Dof.vy], i] = -1; # vy = 0
  end
  for i in bdry_nodes[RIGHT]  # right edge
    ID[dofs[Dof.vy], i] = -1; # vy = 0
  end

  # Neumann boundary conditions
  inh_neu_bcs = [(LEFT, STRETCH, 4.0), (RIGHT, STRETCH, 4.0)];

  return dofs, ndf, ID, inh_dir_bcs, inh_neu_bcs;
end # get_f_coue_bc_info


"""
    get_f_pull_bc_info(numnp, IX, bdry_nodes, bdry_inner_nodes, p::Params; args...)

Generate [`Mesh`](@ref) data for the flat tube-pulling [`Scenario`](@ref)
`F_PULL`.

Relevant unknowns are determined via [`get_dofs`](@ref MembraneAleFem.get_dofs),
depending on the choice of mesh [`Motion`](@ref).
The following boundary conditions are prescribed for all mesh motions:
- ``v_z = 0`` on all boundaries, and all inner boundaries
- ``\\bm{\\bar{F}} = (k_{\\text{b}} / 4)\\bm{\\nu}`` on all boundaries
- ``v_x = 0``, ``v_y = 0`` at center of each edge
- ``v_x = 0``, ``v_y = 0``, ``v_z = \\,```args[:pull_speed]```\\bm{e}_z`` on
  nodes of the central element

The following are prescribed when the [`Motion`](@ref) is not Lagrangian:
- ``v^{\\text{m}}_x = 0``, ``v^{\\text{m}}_y = 0`` at center of each edge
- ``v^{\\text{m}}_x = 0``, ``v^{\\text{m}}_y = 0``,
  ``v^{\\text{m}}_z = \\,```args[:pull_speed}```\\bm{e}_z`` on nodes of the
  central element

The following are prescribed when the [`Motion`](@ref) is either ALE-viscous or
ALE-viscous-bending:
- ``\\bm{v}^{\\text{m}} = \\bm{0}`` on all boundaries

The following are prescribed when the [`Motion`](@ref) is ALE-viscous-bending:
- ``v^{\\text{m}}_z = 0`` on all inner boundaries
"""
function get_f_pull_bc_info(
    numnp::Int64,
    IX::Matrix{Int64},
    bdry_nodes::Dict{Boundary, Vector{Int64}},
    bdry_inner_nodes::Dict{Boundary, Vector{Int64}},
    p::Params;
    args...
  )

  # fundamental unknowns
  dofs = get_dofs(p.motion);
  ndf  = length(values(dofs));

  ## construct ID matrix
  ID = zeros(Int64, ndf, numnp);

  ## remove Dirichlet boundary conditions
  inh_dir_bcs = Tuple{Dof.Unknown, Int64, Float64}[];

  # vz = 0 on all boundaries
  for bdry ∈ keys(bdry_nodes), i ∈ bdry_nodes[bdry]
    ID[dofs[Dof.vz], i] = -1; # vz = 0
    if p.motion == ALEV || p.motion == ALEVB
      # vm = 0 on all boundaries
      ID[dofs[Dof.vmx], i] = -1; # vmx = 0
      ID[dofs[Dof.vmy], i] = -1; # vmy = 0
      ID[dofs[Dof.vmz], i] = -1; # vmz = 0
    end
  end

  # vz = 0 on all inner boundaries
  for bdry ∈ keys(bdry_inner_nodes), i ∈ bdry_inner_nodes[bdry]
    ID[dofs[Dof.vz], i] = -1; # vz = 0
    if p.motion == ALEVB
      # zero slope condition on mesh for ALE-viscous-bending
      ID[dofs[Dof.vmz], i] = -1; # vmz = 0
    end
  end

  # pull tube from the center of the domain
  center_elem_id = get_pull_el_id(size(IX)[2]);
  for pull_node_id in IX[:,center_elem_id]

    ID[dofs[Dof.vx], pull_node_id] = -1;
    ID[dofs[Dof.vy], pull_node_id] = -1;
    ID[dofs[Dof.vz], pull_node_id] = -1;
    push!(inh_dir_bcs, (Dof.vz, pull_node_id, args[:pull_speed]));

    if p.motion != LAG
      ID[dofs[Dof.vmx], pull_node_id] = -1;
      ID[dofs[Dof.vmy], pull_node_id] = -1;
      ID[dofs[Dof.vmz], pull_node_id] = -1;
      push!(inh_dir_bcs, (Dof.vmz, pull_node_id, args[:pull_speed]));
    end

  end


  # prevent rigid body rotations
  for bdry ∈ keys(bdry_nodes)
    bdry_center_node_id = bdry_nodes[bdry][floor(Int64, end/2)+1];
    ID[dofs[Dof.vx], bdry_center_node_id] = -1;
    ID[dofs[Dof.vy], bdry_center_node_id] = -1;

    if p.motion != LAG
      ID[dofs[Dof.vmx], bdry_center_node_id] = -1;
      ID[dofs[Dof.vmy], bdry_center_node_id] = -1;
    end
  end

  # Neumann boundary conditions
  λval  = p.kb / 4;
  inh_neu_bcs = [(LEFT,   STRETCH, λval),
                 (RIGHT,  STRETCH, λval),
                 (TOP,    STRETCH, λval),
                 (BOTTOM, STRETCH, λval)
                ];

  @assert p.pn == 0.0 "$(F_PULL) with a normal pressure is not implemented";

  return dofs, ndf, ID, inh_dir_bcs, inh_neu_bcs;
end # get_f_pull_bc_info


"""
    get_f_pois_bc_info(numnp, bdry_nodes)

Generate [`Mesh`](@ref) data for the flat Hagen--Poiseuille [`Scenario`](@ref)
`F_POIS`.

In this scenario, a `STATIC` mesh [`Motion`](@ref) is required; relevant
[`Unknown`](@ref Dof.Unknown)s are `vx`, `vy`, and `λ`.
The following boundary conditions are prescribed:
- ``v_x = 0.0`` along the top and bottom edges (corners included)
- ``v_y = 0.0`` on all edges
- ``\\bm{\\bar{F}} = 4 \\bm{\\nu}`` on the left edge
- ``\\bm{\\bar{F}} = 8 \\bm{\\nu}`` on the right edge
"""
function get_f_pois_bc_info(
    numnp::Int64,
    bdry_nodes::Dict{Boundary, Vector{Int64}}
  )

  # fundamental unknowns
  dofs = Dict(Dof.vx=>1, Dof.vy=>2, Dof.λ=>3);
  ndf = length(values(dofs));

  ## construct ID matrix
  ID = zeros(Int64, ndf, numnp);

  ## remove Dirichlet boundary conditions
  inh_dir_bcs = Tuple{Dof.Unknown, Int64, Float64}[];

  # top and bottom edges
  for bdry ∈ (TOP, BOTTOM), i ∈ bdry_nodes[bdry]
    ID[dofs[Dof.vx], i] = -1; # vx = 0
    ID[dofs[Dof.vy], i] = -1; # vy = 0
  end

  # left and right edges
  for bdry ∈ (LEFT, RIGHT), i ∈ bdry_nodes[bdry]
    ID[dofs[Dof.vy], i] = -1; # vy = 0
  end

  # Neumann boundary conditions
  inh_neu_bcs = [(LEFT, STRETCH, 4.0), (RIGHT, STRETCH, 8.0)];

  return dofs, ndf, ID, inh_dir_bcs, inh_neu_bcs;
end # get_f_pois_bc_info


"""
    get_f_bend_bc_info(numnp, bdry_nodes, p::Params; args...)

Generate [`Mesh`](@ref) data for the flat bending [`Scenario`](@ref) `F_BEND`.

Relevant unknowns are determined via [`get_dofs`](@ref MembraneAleFem.get_dofs),
depending on the choice of mesh [`Motion`](@ref).
The following boundary conditions are prescribed for all mesh motions:
- `LEFT` edge:
  - ``\\bm{v} = \\bm{0}``
  - ``M = \\,```args[:bend_mf]` (see `Params.jl`)
  - if ALE mesh motion: ``\\bm{v}^{\\text{m}} = \\bm{0}``
- `RIGHT` edge:
  - ``f_x = 0`` and ``f_y = 0``
  - ``v_z = 0``
  - ``M = \\,```args[:bend_mf]` (see `Params.jl`)
  - if ALE mesh motion: ``v^{\\text{m}}_z = 0``
- `BOTTOM` edge:
  - ``f_x = 0`` and ``f_z = 0``
  - ``v_y = 0``
  - ``M = 0``
  - if ALE mesh motion: ``v^{\\text{m}}_y = 0``
- `TOP` edge:
  - ``f_x = 0`` and ``f_z = 0``
  - ``v_y = 0``
  - ``M = 0``
  - if ALE mesh motion: ``v^{\\text{m}}_y = 0``
"""
function get_f_bend_bc_info(
    numnp::Int64,
    bdry_nodes::Dict{Boundary, Vector{Int64}},
    p::Params;
    args...
  )

  # fundamental unknowns
  dofs = get_dofs(p.motion);
  ndf  = length(values(dofs));

  ## construct ID matrix
  ID = zeros(Int64, ndf, numnp);

  ## remove Dirichlet boundary conditions
  inh_dir_bcs = Tuple{Dof.Unknown, Int64, Float64}[];

  # top and bottom edges
  for bdry ∈ (TOP, BOTTOM), i ∈ bdry_nodes[bdry] # bottom edge
    ID[dofs[Dof.vy], i] = -1; # vy = 0
    if p.motion == ALEV || p.motion == ALEVB
      ID[dofs[Dof.vmy], i] = -1; # vmy = 0
    end
  end
  for i ∈ bdry_nodes[LEFT]   # left edge
    ID[dofs[Dof.vx], i] = -1; # vx = 0
    ID[dofs[Dof.vy], i] = -1; # vy = 0
    ID[dofs[Dof.vz], i] = -1; # vz = 0
    if p.motion == ALEV || p.motion == ALEVB
      ID[dofs[Dof.vmx], i] = -1; # vmx = 0
      ID[dofs[Dof.vmy], i] = -1; # vmy = 0
      ID[dofs[Dof.vmz], i] = -1; # vmz = 0
    end
  end
  for i ∈ bdry_nodes[RIGHT]  # right edge
    ID[dofs[Dof.vz], i] = -1; # vz = 0
    if p.motion == ALEV || p.motion == ALEVB
      ID[dofs[Dof.vmz], i] = -1; # vmz = 0
    end
  end

  # Neumann boundary conditions;
  inh_neu_bcs = [(LEFT,  MOMENT, args[:bend_mf]),
                 (RIGHT, MOMENT, args[:bend_mf])];

  return dofs, ndf, ID, inh_dir_bcs, inh_neu_bcs;
end # get_f_bend_bc_info


"""
    get_dofs(motion::Motion)::Dict{Dof.Unknown, Int64}

Return global ordering of [`Unknown`](@ref Dof.Unknown) degrees of freedom,
depending on the mesh [`Motion`](@ref).

The following is the mapping from motion to unknowns:
- Lagrangian --> ``v_x``, ``v_y``, ``v_z``, ``λ``
- Eulerian --> ``v_x``, ``v_y``, ``v_z``, ``v^{\\text{m}}_x``,
  ``v^{\\text{m}}_y``, ``v^{\\text{m}}_z``, ``λ``
- ALE --> ``v_x``, ``v_y``, ``v_z``, ``v^{\\text{m}}_x``, ``v^{\\text{m}}_y``,
  ``v^{\\text{m}}_z``, ``λ``, ``p^{\\text{m}}``

Note that the flat, 2-D [`Scenario`](@ref)s (`F_CAVI`, `F_COUE`, `F_POIS`) have
a reduced number of degrees of freedom, because there are no unknowns in the
``z``-direction.
Thus, this function is not called in those scenarios.
"""
function get_dofs(
    motion::Motion
  )::Dict{Dof.Unknown, Int64}

  if motion == LAG || motion == STATIC
    return Dict(Dof.vx=>1,   Dof.vy=>2,  Dof.vz=>3,  Dof.λ=>4);
  elseif motion == EUL
    return Dict(Dof.vx=>1,   Dof.vy=>2,  Dof.vz=>3,
                Dof.vmx=>4,  Dof.vmy=>5, Dof.vmz=>6, Dof.λ=>7);
  elseif motion == ALEV || motion == ALEVB
    return Dict(Dof.vx=>1,   Dof.vy=>2,  Dof.vz=>3,
                Dof.vmx=>4,  Dof.vmy=>5, Dof.vmz=>6, Dof.λ=>7, Dof.pm=>8);
  else
    @assert false "$(motion) motion degrees of freedom not provided";
    return Nothing;
  end

end # get_dofs

