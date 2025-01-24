

export Mesh
export get_basis_fns, get_N


"""
    Mesh(p::Params; args...)

The generated mesh, which includes the degree-of-freedom numbering and all basis
functions for the given [`Scenario`](@ref).

This struct contains the following scenario-independent fields:
- `num1el` --> number of elements in ``ζ^1`` direction
- `num2el` --> number of elements in ``ζ^2`` direction
- `numel` --> number of area elements
- `num1np` --> number of nodal points in ``ζ^1`` direction
- `num2np` --> number of nodal points in ``ζ^2`` direction
- `numnp` -->  number of nodal points on the 2-d mesh
- `IX` --> matrix containing nonzero node numbers for each area element
- `bdry_elems` --> mapping from [`Boundary`](@ref) to element numbers
- `crnr_elems` --> mapping from [`Corner`](@ref) to element number
- `bdry_nodes` --> mapping from [`Boundary`](@ref) to node numbers
- `bdry_inner_nodes` --> mapping from [`Boundary`](@ref) to inner node numbers
- `crnr_nodes` --> mapping from [`Corner`](@ref) to node number
- `crnr_inner_nodes` --> mapping from [`Corner`](@ref) to inner node number
- `kv1` --> [`KnotVector`](@ref) along ``ζ^1`` direction
- `kv2` --> [`KnotVector`](@ref) along ``ζ^2`` direction
- `area_gp_fns` --> [`AreaGpBasisFns`](@ref): 2-D area basis functions
- `bdry_gp_fns` --> mapping from [`Boundary`](@ref) to [`BdryGpBasisFns`](@ref)
- `crnr_gp_fns` --> mapping from [`Corner`](@ref) to 2-D basis functions
  [`GpBasisFnsζα`](@ref)

The generated mesh depends on the problem being solved.
The function [`generate_scenario`](@ref MembraneAleFem.generate_scenario) uses
the scenario-specific results from `Bc.jl` to set the following
scenario-dependent fields:
- `topology` --> surface [`Topology`](@ref)
- `dofs` --> list of active [`Unknown`](@ref Dof.Unknown)s with global ordering
- `ndf` --> number of active [`Unknown`](@ref Dof.Unknown)s
- `ID` --> matrix mapping node number to indices of unknowns at that node
- `ID_inv` --> inverse of `ID`: map `unknown_id` to `(node_number, dof_number)`
- `nmdf` --> total number of mesh degrees of freedom
- `LM` --> map from area element to indices of all unknowns of that element
- `inh_dir_bcs` --> list of inhomogeneous Dirichlet boundary conditions
- `inh_neu_bcs` --> list of inhomogeneous Neumann boundary conditions
"""
struct Mesh
  # scenario independent
  num1el::Int64
  num2el::Int64
  numel::Int64
  num1np::Int64
  num2np::Int64
  numnp::Int64
  IX::Matrix{Int64}
  bdry_elems::Dict{Boundary, Vector{Int64}}
  crnr_elems::Dict{Corner, Int64}
  bdry_nodes::Dict{Boundary, Vector{Int64}}
  bdry_inner_nodes::Dict{Boundary, Vector{Int64}}
  crnr_nodes::Dict{Corner, Int64}
  crnr_inner_nodes::Dict{Corner, Int64}
  kv1::KnotVector
  kv2::KnotVector
  area_gp_fns::AreaGpBasisFns
  bdry_gp_fns::Dict{Boundary, BdryGpBasisFns}
  crnr_gp_fns::Dict{Corner, GpBasisFnsζα}
  # scenario dependent
  topology::Topology
  dofs::Dict{Dof.Unknown, Int64}
  ndf::Int64
  ID::Matrix{Int64}
  ID_inv::Vector{CartesianIndex{2}}
  nmdf::Int64
  LM::Matrix{Int64}
  inh_dir_bcs::Vector{Tuple{Dof.Unknown, Int64, Float64}}
  inh_neu_bcs::Vector{Tuple{Boundary, Neumann, Float64}}

  Mesh(p::Params; args...) = new(generate_mesh(p; args...)...);
end # Mesh


"""
    generate_mesh(p::Params; args...)

External constructor of the `Mesh` `struct`.

This function carries out all scenario-independent calculations, and then
calls the appropriate function for each scenario.

See also
[`Mesh`](@ref)
"""
function generate_mesh(
    p::Params;
    args...
  )

  if p.output
    io = open(args[:out_path] * "/" * args[:out_file], "a");
    println(io, "-> creating a ", p.num1el, " x ", p.num2el,
                " mesh for the ", p.scenario, " scenario");
  end

  topology = get_topology(p.scenario);

  numαel = (p.num1el, p.num2el);
  numel =  prod(numαel);         # number of area elements

  if topology == FLAT
    num1np = p.num1el + p.poly;  # nodes in ζ1 = x/p.length direction
    num2np = p.num2el + p.poly;  # nodes in ζ2 = y/p.length direction

  elseif topology == CYLINDER
    num1np = p.num1el;           # nodes in ζ1 = θ/2π       direction
    num2np = p.num2el + p.poly;  # nodes in ζ2 = z/p.length direction

  else
    @assert false "mesh construction for $(topology) topology not implemented";
  end

  numαnp = (num1np, num2np);
  numnp  = prod(numαnp);     # total number of nodal points

  # boundary element indices on each boundary
  bdry_elems = Dict{Boundary, Vector{Int64}}(
    BOTTOM => [i for i=1:p.num1el],
    RIGHT  => [i for i=p.num1el:p.num1el:numel],
    TOP    => [i for i=numel-p.num1el+1:numel],
    LEFT   => [i for i=1:p.num1el:numel-p.num1el+1]
  );

  # element indices of each corner
  crnr_elems = Dict{Corner, Int64}(
    BOTTOM_LEFT  => 1,
    BOTTOM_RIGHT => p.num1el,
    TOP_LEFT     => numel-p.num1el+1,
    TOP_RIGHT    => numel
  );

  if topology == FLAT

    # nodal indices on each boundary
    bdry_nodes = Dict{Boundary, Vector{Int64}}(
      BOTTOM => [i for i=1:num1np],
      RIGHT  => [i for i=num1np:num1np:numnp],
      TOP    => [i for i=numnp-num1np+1:numnp],
      LEFT   => [i for i=1:num1np:numnp-num1np+1]
    );

    # inner nodal indices on each boundary
    bdry_inner_nodes = Dict{Boundary, Vector{Int64}}(
      BOTTOM => [i+num1np for i=1:num1np],
      RIGHT  => [i-1      for i=num1np:num1np:numnp],
      TOP    => [i-num1np for i=numnp-num1np+1:numnp],
      LEFT   => [i+1      for i=1:num1np:numnp-num1np+1]
    );

    # nodal indices of each corner
    crnr_nodes = Dict{Corner, Int64}(
      BOTTOM_LEFT  => 1,
      BOTTOM_RIGHT => num1np,
      TOP_LEFT     => numnp-num1np+1,
      TOP_RIGHT    => numnp
    );

    # inner nodal indices of each corner
    crnr_inner_nodes = Dict{Corner, Int64}(
      BOTTOM_LEFT  => num1np+2,
      BOTTOM_RIGHT => 2*num1np-1,
      TOP_LEFT     => numnp-2*num1np+2,
      TOP_RIGHT    => numnp-num1np-1
    );

    # boundary knot vectors
    if p.scenario == F_PULL && p.num1el >= 18 && p.num2el >= 18
      if p.output
        println(io, "constructing a fine mesh at the center");
      end
      kv1 = KnotVector(get_fine_ζs(p.num1el, p.poly), p.poly, CLAMPED);
      kv2 = KnotVector(get_fine_ζs(p.num2el, p.poly), p.poly, CLAMPED);
    else
      kv1 = KnotVector(p.num1el, p.poly, CLAMPED);
      kv2 = KnotVector(p.num2el, p.poly, CLAMPED);
    end

  elseif topology == CYLINDER

    # nodal indices on each boundary
    bdry_nodes = Dict{Boundary, Vector{Int64}}(
      BOTTOM => [i for i=1:num1np],
      TOP    => [i for i=numnp-num1np+1:numnp],
    );
    # inner nodal indices on each boundary
    bdry_inner_nodes = Dict{Boundary, Vector{Int64}}(
      BOTTOM => [i+num1np for i=1:num1np],
      TOP    => [i-num1np for i=numnp-num1np+1:numnp],
    );

    # no corner nodes in a cylindrical mesh
    crnr_nodes = Dict{Corner, Int64}();
    crnr_inner_nodes = Dict{Corner, Int64}();

    # boundary knot vectors
    kv1 = KnotVector(p.num1el, p.poly, CLOSED);  # periodic in θ
    kv2 = KnotVector(p.num2el, p.poly, CLAMPED);

  else
    @assert false "mesh construction for $(topology) topology not implemented";
  end

  # 1-D line basis functions
  line_gp_fns1 = LineGpBasisFns(kv1, p.gp1d);
  line_gp_fns2 = LineGpBasisFns(kv2, p.gp1d);

  # area basis functions
  area_gp_fns  = AreaGpBasisFns(line_gp_fns1, line_gp_fns2);

  # 2-D boundary basis functions
  bdry_gp_fns = Dict{Boundary, BdryGpBasisFns}(
    BOTTOM => BdryGpBasisFns(line_gp_fns1, line_gp_fns2.ζmin_fns, BOTTOM),
    RIGHT  => BdryGpBasisFns(line_gp_fns2, line_gp_fns1.ζmax_fns, RIGHT),
    TOP    => BdryGpBasisFns(line_gp_fns1, line_gp_fns2.ζmax_fns, TOP),
    LEFT   => BdryGpBasisFns(line_gp_fns2, line_gp_fns1.ζmin_fns, LEFT)
  );

  # corner basis functions
  crnr_gp_fns = Dict{Corner, GpBasisFnsζα}(
    BOTTOM_LEFT  => GpBasisFnsζα(line_gp_fns1.ζmin_fns, line_gp_fns2.ζmin_fns),
    BOTTOM_RIGHT => GpBasisFnsζα(line_gp_fns1.ζmax_fns, line_gp_fns2.ζmin_fns),
    TOP_LEFT     => GpBasisFnsζα(line_gp_fns1.ζmin_fns, line_gp_fns2.ζmax_fns),
    TOP_RIGHT    => GpBasisFnsζα(line_gp_fns1.ζmax_fns, line_gp_fns2.ζmax_fns)
  );

  # construct IX matrix
  IX = construct_IX(kv1, kv2, num1np);

  if p.output
    close(io);
  end

  # messy way to return all relevant information
  return p.num1el, p.num2el, numel, num1np, num2np, numnp, IX,
    bdry_elems, crnr_elems, bdry_nodes, bdry_inner_nodes, crnr_nodes,
    crnr_inner_nodes, kv1, kv2, area_gp_fns, bdry_gp_fns, crnr_gp_fns,
    topology, generate_scenario(numel, numnp, IX, bdry_nodes,
                                bdry_inner_nodes, crnr_nodes, p; args...)...;
end # generate_mesh


"""
    generate_scenario(numel, numnp, IX, bdry_nodes, ..., p::Params; args...)

Apply scenario-specific boundary conditions towards mesh generation.

Checks to see whether a helper function is available in `Bc.jl` for the given
[`Scenario`](@ref), and then constructs the `ID` matrix that maps from nodes to
global degrees of freedom.
The inverse `ID_inv` is also generated, as is the `LM` matrix—which provides the
local degrees of freedom for each element.
"""
function generate_scenario(
    numel::Int64,
    numnp::Int64,
    IX::Matrix{Int64},
    bdry_nodes::Dict{Boundary, Vector{Int64}},
    bdry_inner_nodes::Dict{Boundary, Vector{Int64}},
    crnr_nodes::Dict{Corner, Int64},
    p::Params;
    args...
  )

  dofs, ndf, ID, inh_dir_bcs, inh_neu_bcs = get_scenario_bc_info(numnp, IX,
                                              bdry_nodes, bdry_inner_nodes,
                                              crnr_nodes, p; args...);
  dof_index = 1;
  for node_id = 1:numnp, dof_id = 1:ndf
    if ID[dof_id, node_id] != -1
      ID[dof_id, node_id] = dof_index;
      dof_index += 1;
    else
      ID[dof_id, node_id] = 0;
    end
  end

  # total number of degrees of freedom, after Dirichlet boundary conditions have
  # been removed
  nmdf = maximum(ID);

  # create inverse ID mapping: unknown_id => node_id, dof_id
  ID_inv = CartesianIndex{2}[];
  for node_id = 1:numnp, dof_id=1:ndf
    if ID[dof_id, node_id] != 0
      push!(ID_inv, CartesianIndex(node_id, dof_id));
    end
  end

  # construct LM matrices
  LM = reshape(ID[:,IX], (ndf*NEN, numel));

  return dofs, ndf, ID, ID_inv, nmdf, LM, inh_dir_bcs, inh_neu_bcs;
end # generate_scenario


"""
    get_basis_fns(el_id::Int64, gp_id::Int64, mesh::Mesh)::GpBasisFnsζα

Return [`GpBasisFnsζα`](@ref) of the `gp_id```^{\\text{th}}`` 2-D Gauss point of
the `el_id```^{\\text{th}}`` area element.
"""
function get_basis_fns(
    el_id::Int64,    # area element id
    gp_id::Int64,    # 2-D
    mesh::Mesh
  )::GpBasisFnsζα

  @assert (gp_id >= 1 && gp_id <= GP1D^2) "2-D Gauss point index out of bounds";
  return mesh.area_gp_fns.ufns[mesh.area_gp_fns.uel_ids[el_id], gp_id];
end # get_basis_fns


"""
    get_gpw(el_id::Int64, gp_id::Int64, mesh::Mesh)::Float64

Return the Gauss point weight of the call to [`get_basis_fns`](@ref).
"""
function get_gpw(
    el_id::Int64,    # area element id
    gp_id::Int64,    # 2-D
    mesh::Mesh
  )::Float64

  return get_basis_fns(el_id, gp_id, mesh).w;
end # get_gpw


"""
    get_N(el_id::Int64, gp_id::Int64, mesh::Mesh)::SVector{NEN,Float64}

Return local basis functions ``[𝐍^e]`` of the call to [`get_basis_fns`](@ref).
"""
function get_N(
    el_id::Int64,    # area element id
    gp_id::Int64,    # 2-D
    mesh::Mesh
  )::SVector{NEN,Float64}

  return get_basis_fns(el_id, gp_id, mesh).N;
end # get_N


"""
    get_∂Nα(el_id::Int64, gp_id::Int64, mesh::Mesh)::SMatrix{NEN,ζDIM,Float64}

Return local basis function derivatives ``[𝐍^e]_{, α}`` of the call to
[`get_basis_fns`](@ref).
"""
function get_∂Nα(
    el_id::Int64,    # area element id
    gp_id::Int64,    # 2-D
    mesh::Mesh
  )::SMatrix{NEN,ζDIM,Float64}

  return get_basis_fns(el_id, gp_id, mesh).∂Nα;
end # get_∂Nα


"""
    get_∂∂Nαβ(el_id::Int64, gp_id::Int64, mesh::Mesh)::SMatrix{NEN,VOIGT,Float64}

Return local basis function derivatives ``[𝐍^e]_{, α β}`` of the call to
[`get_basis_fns`](@ref).
"""
function get_∂∂Nαβ(
    el_id::Int64,    # area element id
    gp_id::Int64,    # 2-D
    mesh::Mesh
  )::SMatrix{NEN,VOIGT,Float64}

  return get_basis_fns(el_id, gp_id, mesh).∂∂Nαβ;
end # get_∂∂Nαβ


"""
    get_basis_fns(bdry::Boundary, el_id::Int64, gp_id::Int64, mesh::Mesh)::GpBasisFnsζα

Return [`GpBasisFnsζα`](@ref) of the `gp_id```^{\\text{th}}`` 1-D Gauss point of
the `el_id```^{\\text{th}}`` element on the [`Boundary`](@ref) `bdry`.
"""
function get_basis_fns(
    bdry::Boundary,
    el_id::Int64,    # area element id that is on the boundary
    gp_id::Int64,    # 1-D
    mesh::Mesh
  )::GpBasisFnsζα

  @assert (gp_id >= 1 && gp_id <= GP1D) "1-D Gauss point index out of bounds";
  bel_id = findfirst(==(el_id), mesh.bdry_elems[bdry]);
  @assert !isequal(typeof(bel_id), Nothing) "element id not found on boundary";

  return mesh.bdry_gp_fns[bdry].ufns[
           mesh.bdry_gp_fns[bdry].uel_ids[bel_id], gp_id];
end # get_basis_fns


"""
    get_gpw(bdry::Boundary, el_id::Int64, gp_id::Int64, mesh::Mesh)::Vector{Float64}

Return the Gauss point weight of the call to [`get_basis_fns`](@ref).
"""
function get_gpw(
    bdry::Boundary,
    el_id::Int64,    # area element id that is on the boundary
    gp_id::Int64,    # 1-D
    mesh::Mesh
  )::Float64

  return get_basis_fns(bdry, el_id, gp_id, mesh).w;
end # get_gpw


"""
    get_N(bdry::Boundary, el_id::Int64, gp_id::Int64, mesh::Mesh)

Return local basis functions ``[𝐍^e]`` of the call to [`get_basis_fns`](@ref).
"""
function get_N(
    bdry::Boundary,
    el_id::Int64,    # area element id that is on the boundary
    gp_id::Int64,    # 1-D
    mesh::Mesh
  )::SVector{NEN,Float64}

  return get_basis_fns(bdry, el_id, gp_id, mesh).N;
end # get_N


"""
    get_∂Nα(bdry::Boundary, el_id::Int64, gp_id::Int64, mesh::Mesh)

Return local basis function derivatives ``[𝐍^e]_{, α}`` of the call to
[`get_basis_fns`](@ref).
"""
function get_∂Nα(
    bdry::Boundary,
    el_id::Int64,    # area element id that is on the boundary
    gp_id::Int64,    # 1-D
    mesh::Mesh
  )::SMatrix{NEN,ζDIM,Float64}

  return get_basis_fns(bdry, el_id, gp_id, mesh).∂Nα;
end # get_∂Nα


"""
    get_∂∂Nαβ(bdry::Boundary, el_id::Int64, gp_id::Int64, mesh::Mesh)

Return local basis function derivatives ``[𝐍^e]_{, α β}`` of the call to
[`get_basis_fns`](@ref).
"""
function get_∂∂Nαβ(
    bdry::Boundary,
    el_id::Int64,    # area element id that is on the boundary
    gp_id::Int64,    # 1-D
    mesh::Mesh
  )::SMatrix{NEN,VOIGT,Float64}

  return get_basis_fns(bdry, el_id, gp_id, mesh).∂∂Nαβ;
end # get_∂∂Nαβ


"""
    get_v_order(Mesh.dofs)::Vector{Int64}

Return active velocity degree of freedom indices.
"""
function get_v_order(
    dofs::Dict{Dof.Unknown, Int64}
  )::Vector{Int64}
  return [get(dofs, Dof.vx, 0);
          get(dofs, Dof.vy, 0);
          get(dofs, Dof.vz, 0)];
end # get_v_order


"""
    get_m_order(Mesh.dofs)::Vector{Int64}

Return active mesh velocity degree of freedom indices.
"""
function get_m_order(
    dofs::Dict{Dof.Unknown, Int64}
  )::Vector{Int64}
  return [get(dofs, Dof.vmx, 0);
          get(dofs, Dof.vmy, 0);
          get(dofs, Dof.vmz, 0)];
end # get_m_order


"""
    get_λ_order(Mesh.dofs)::Int64

Return surface tension degree of freedom index.
"""
function get_λ_order(
    dofs::Dict{Dof.Unknown, Int64}
  )::Int64
  return get(dofs, Dof.λ, 0);
end # get_λ_order


"""
    get_p_order(Mesh.dofs)::Int64

Return mesh pressure degree of freedom index.
"""
function get_p_order(
    dofs::Dict{Dof.Unknown, Int64}
  )::Int64
  return get(dofs, Dof.pm, 0);
end # get_p_order


"""
    get_m_motion_order(motion::Motion, Mesh.dofs)::Vector{Int64}

Return degree of freedom indices corresponding to the mesh motion.
"""
function get_m_motion_order(
    motion::Motion,
    dofs::Dict{Dof.Unknown, Int64}
  )::Vector{Int64}

  if motion == STATIC
    return [0; 0; 0];
  elseif motion == LAG
    return get_v_order(dofs);
  else
    return get_m_order(dofs);
  end

end # get_m_motion_order


"""
    get_topology(scenario::Scenario)::Topology

Return mesh `topology` corresponding to a given `scenario`.
"""
function get_topology(
    scenario::Scenario
  )::Topology

  flat_scenarios = [F_CAVI, F_COUE, F_POIS, F_PULL, F_BEND];
  cyl_scenarios  = [];

  if scenario ∈ flat_scenarios
    return FLAT;
  elseif scenario ∈ cyl_scenarios
    return CYLINDER;
  else
    @assert false "topology for $(scenario) scenario not specified";
    return Nothing;
  end

end


"""
construct_IX()

Construct the `IX` matrix.
"""
function construct_IX(
    kv1::KnotVector,
    kv2::KnotVector,
    num1np::Int64
  )::Matrix{Int64}

  IX = zeros(Int64, NEN, kv1.nel * kv2.nel);

  for eid2 = 1:kv2.nel, eid1 = 1:kv1.nel
    spline_ids2 = get_bspline_indices(kv2, eid2+POLY);
    spline_ids1 = get_bspline_indices(kv1, eid1+POLY);

    for nid2 = 1:POLY+1, nid1 = 1:POLY+1
      IX[nid1 + (nid2-1)*(POLY+1), eid1 + (eid2-1)*kv1.nel] =
        spline_ids1[nid1] + num1np * (spline_ids2[nid2] - 1);
    end
  end

  return IX;
end # construct_IX



