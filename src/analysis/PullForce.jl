

"""
    get_pull_el_id(numel::Int64)

Return the id of the element that will be pulled, given the number of elements.
"""
function get_pull_el_id(
    numel::Int64
  )

  return ceil(Int64, numel / 2);
end # get_pull_el_id


"""
    get_adj_maps(num1el::Int64, numel::Int64, IX::Matrix{Int64})

For each element adjacent to the pulled one, return local ids of pulled nodes.

These mappings are required to determine the residual vector of pulled nodes,
even though the nodes are treated as Dirichlet boundary conditions and thus
excluded from the list of degrees of freedom.
"""
function get_adj_maps(
    num1el::Int64,
    numel::Int64,
    IX::Matrix{Int64},
    poly::Int64
  )

  # element over which pull force is applied/calculated
  pull_el_id = get_pull_el_id(numel);
  pull_nodes = IX[:, pull_el_id];

  # adjacent elements contributing to pull force
  adj_el_ids = vcat([[pull_el_id+i*num1el-poly:pull_el_id+i*num1el+poly;]
                        for i=-poly:poly]...);

  # local node indices within local residual vectors (adjacent, pull)
  adj_node_map = Tuple{Vector{Int64}, Vector{Int64}}[];
  for adj_el_id in adj_el_ids
    adj_nodes = IX[:, adj_el_id];
    push!(adj_node_map, (
      vcat( [(i-1)*XDIM .+ [1:XDIM;]
             for i in findall(x->x in adj_nodes, pull_nodes)]...),
      vcat( [(i-1)*XDIM .+ [1:XDIM;]
             for i in findall(x->x in pull_nodes, adj_nodes)]...))
    );
  end

  return adj_el_ids, adj_node_map;
end # get_adj_maps


"""
    calc_pull_force(mesh, xms, cps, adj_el_ids, adj_node_map, p::Params)

Return the calculated pull force 𝓕, as a 3×1 vector.
"""
function calc_pull_force(
    mesh::Mesh,
    xms::Matrix{Float64},
    cps::Matrix{Float64},
    adj_el_ids::Vector{Int64},
    adj_node_map::Vector{Tuple{Vector{Int64}, Vector{Int64}}},
    p::Params
  )::Vector{Float64}

  rv_pull = zeros(Float64, XDIM*NEN);

  for (id, adj_el_id) in enumerate(adj_el_ids)
    xms_el = complex(xms[mesh.IX[:,adj_el_id], :]);
    cps_el = complex(cps[mesh.IX[:,adj_el_id], :]);
    rv_el, _, _, _ = calc_elem_dof_residuals(mesh, adj_el_id, xms_el, cps_el, p);
    rv_pull[adj_node_map[id][1]] += real(rv_el[adj_node_map[id][2]]);
  end

  return vec(sum(reshape(rv_pull, XDIM, :), dims=2));
end


