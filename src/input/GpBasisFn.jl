

export GpBasisFnsζ, GpBasisFnsζα, LineGpBasisFns, BdryGpBasisFns, AreaGpBasisFns


"""
    GpBasisFnsζ(gpwζ::Float64, ζ::Float64, kv::KnotVector)

1-D basis functions ``N(ζ)`` at the provided quadrature point `ζ` with weight
`gpwζ`.

The B-spline values and derivatives are determined with
[`get_bspline_ders`](@ref) and stored here.
The struct contains four fields:
- `w`   --> Gauss point weight, passed in as `gpw`
- `N`   --> (`poly`+1)×1 vector of nonzero basis functions ``N(ζ)``
- `dN`  --> (`poly`+1)×1 vector of nonzero first derivatives ``N'(ζ)``
- `ddN` --> (`poly`+1)×1 vector of nonzero second derivatives ``N''(ζ)``
"""
struct GpBasisFnsζ
  w::Float64
  N::SVector{POLY+1,Float64}
  dN::SVector{POLY+1,Float64}
  ddN::SVector{POLY+1,Float64}

  GpBasisFnsζ(
    gpwζ::Float64,
    ζ::Float64,
    kv::KnotVector
  ) = new(
    gp_basis_fns_ζ(gpwζ, ζ, kv)...
  );

end # GpBasisFnsζ


"""
    gp_basis_fns_ζ(gpw::Float64, ζ::Float64, kv::KnotVector)

External constructor of the `GpBasisFnsζ` struct.

See also
[`GpBasisFnsζ`](@ref)
"""
function gp_basis_fns_ζ(
    gpwζ::Float64,
    ζ::Float64,
    kv::KnotVector
  )
  bspline_ders = get_bspline_ders(kv, ζ, NDERS);
  return gpwζ,
         SVector{POLY+1,Float64}(bspline_ders[:,1]),
         SVector{POLY+1,Float64}(bspline_ders[:,2]),
         SVector{POLY+1,Float64}(bspline_ders[:,3]);
end


# override == to compare each field
# from https://discourse.julialang.org/t/comparing-julia-structs/101554/2
Base.:(==)(fn1::GpBasisFnsζ, fn2::GpBasisFnsζ) =
  isequal(fn1.w,   fn2.w)  &&
  isequal(fn1.N,   fn2.N)  &&
  isequal(fn1.dN,  fn2.dN) &&
  isequal(fn1.ddN, fn2.ddN);


# hash must also be altered, since == is used to determine if keys are equal
function Base.hash(gpfnζ::GpBasisFnsζ, h_in::UInt)
  ht = hash(GpBasisFnsζ, h_in);
  h1 = hash(gpfnζ.w,   ht);
  h2 = hash(gpfnζ.N,   h1);
  h3 = hash(gpfnζ.dN,  h2);
  h4 = hash(gpfnζ.ddN, h3);
  return h4;
end


"""
    GpBasisFnsζα(fns1::GpBasisFnsζ, fns2::GpBasisFnsζ)

2-D basis functions ``N(ζ^α)`` at a single Gauss point ``ζ^α``, as set by 1-D
[`GpBasisFnsζ`](@ref) `fns1` and `fns2`.

Here ``ζ^1`` and ``ζ^2`` respectively correspond to the `ζ` values passed to
`fns1` and `fns2`, even though this information is not stored.
The 2-D basis functions are calculated as the tensor product of 1-D functions.
The struct contains four fields:
- `w` --> Gauss point weight, given by `fns1.w × fns2.w`
- `N` --> `NEN`×1 vector of nonzero basis functions ``N(ζ^α)``, ordered by the
  tensor product structure
- `∂Nα` --> `NEN`×2 matrix of first derivatives: rows respectively contain
  ``N_{, 1}`` and ``N_{, 2}``
- `∂∂Nαβ` --> `NEN`×3 matrix of second derivatives: rows respectively contain
  ``N_{, 1 1}``, ``N_{, 2 2}``, and ``N_{, 1 2}``
"""
struct GpBasisFnsζα
  w::Float64
  N::SVector{NEN,Float64}
  ∂Nα::SMatrix{NEN,ζDIM,Float64}
  ∂∂Nαβ::SMatrix{NEN,VOIGT,Float64}

  GpBasisFnsζα(
    fns1::GpBasisFnsζ,
    fns2::GpBasisFnsζ,
  ) = new(
    fns1.w * fns2.w,
    SVector{NEN,Float64}([N1*N2 for N2 in fns2.N for N1 in fns1.N]),
    SMatrix{NEN,ζDIM,Float64}([[dN1*N2  for N2  in fns2.N  for dN1 in fns1.dN] [ N1*dN2 for dN2 in fns2.dN for N1  in fns1.N]]),
    SMatrix{NEN,VOIGT,Float64}([[ddN1*N2 for N2 in fns2.N for ddN1 in fns1.ddN] [N1*ddN2 for ddN2 in fns2.ddN for N1 in fns1.N] [dN1*dN2 for dN2 in fns2.dN for dN1 in fns1.dN]])
  );

end # GpBasisFnsζα


"""
    LineGpBasisFns(kv::KnotVector, ngp::Int64)

1-D B-spline basis functions and derivatives at all quadrature points on the
domain of the [`KnotVector`](@ref) `kv`.

Basis functions are determined at each of the `ngp` Gauss points, over each
element.
When using B-splines, scenarios often arise where the basis functions are
identical across many elements.
We accordingly store only the 1-D basis functions over unique elements, as well
as a mapping from elements to unique elements.
The basis functions at the edges of the 1-D parametric domain are stored
separately, as they are required subsequently for application of boundary
conditions.
The struct has five fields:
- `nel` --> number of elements in the [`KnotVector`](@ref) `kv`
- `uel_ids` --> mapping from `elem_id` to `unique_elem_id`
- `ufns` --> `num_unique_elems`×`ngp` matrix of unique 1-D basis functions, of
  type [`GpBasisFnsζ`](@ref)
- `ζmin_fns` --> 1-D basis functions at the min value of the parametric domain
- `ζmax_fns` --> 1-D basis functions at the max value of the parametric domain

Note that this struct contains no information about how basis functions vary in
the orthogonal parametric direction, and thus cannot be directly used to apply
boundary conditions.
See [`BdryGpBasisFns`](@ref).
"""
struct LineGpBasisFns
  nel::Int64
  uel_ids::Vector{Int64}
  ufns::Matrix{GpBasisFnsζ}
  ζmin_fns::GpBasisFnsζ
  ζmax_fns::GpBasisFnsζ

  LineGpBasisFns(
    kv::KnotVector,
    ngp::Int64
  ) = new(
    line_gp_basis_fns(kv, ngp)...
  );
end # LineGpBasisFns


"""
    line_gp_basis_fns(kv::KnotVector, ngp::Int64)

External constructor of the `LineGpBasisFns` struct.

Returns basis functions and `nders` derivatives at each of the `ngp` quadrature
points, for every element along a line (boundary).

See also
[`LineGpBasisFns`](@ref)
"""
function line_gp_basis_fns(
    kv::KnotVector,
    ngp::Int64
  )

  # get unique elements, and the mapping to them from elements
  num_uel, num_el, uel_ids, uel_list = get_unique_1d_elements(kv);

  gpsξ = GaussPointsξ(ngp);
  ugp_fns = Matrix{GpBasisFnsζ}(undef, num_uel, ngp);

  # unique element loop
  for uel_id = 1:num_uel
    gps = GaussPointsζ(gpsξ, uel_list[uel_id][1], uel_list[uel_id][2]);
    # Gauss point loop
    for gp_id = 1:ngp
      ugp_fns[uel_id, gp_id] = GpBasisFnsζ(gps.ws[gp_id], gps.ζs[gp_id], kv);
    end
  end

  # basis functions at ζ_min, ζ_max
  if kv.curve == CLAMPED
    ζmin_fns = GpBasisFnsζ(1.0, kv.ζs[1],   kv);
    ζmax_fns = GpBasisFnsζ(1.0, kv.ζs[end], kv);
  elseif kv.curve == CLOSED
    ζmin_fns = GpBasisFnsζ(1.0, kv.ζs[kv.poly+1],   kv);
    ζmax_fns = GpBasisFnsζ(1.0, kv.ζs[end-kv.poly], kv);
  else
    @assert false "edge basis functions for $(kv.curve) curve not implemented";
  end

  return num_el, uel_ids, ugp_fns, ζmin_fns, ζmax_fns;
end #line_gp_basis_fns


"""
    BdryGpBasisFns(line_gp_fns::LineGpBasisFns, perp_edge_fns::GpBasisFnsζ, bdry::Boundary)

2-D B-spline basis functions and derivatives at all quadrature points on a
boundary.

The [`LineGpBasisFns`](@ref) struct `line_gp_fns` contains all 1-D basis
functions and derivatives along the boundary, with `ngp` Gauss points for each
element.
Here `perp_edge_fns` are the basis functions and derivatives at the boundary in
the orthogonal parametric direction, as contained in a [`GpBasisFnsζ`](@ref)
struct.
For example, on the `RIGHT` boundary, `line_gp_fns` captures 1-D derivatives in
the ``ζ^2`` direction, while `perp_edge_fns` contains 1-D derivatives in the
``ζ^1`` direction on the right edge.
With knowledge of the underlying tensor product structure
[piegl-tiller](@citep), the 2-D basis functions along the specified
[`Boundary`](@ref) `bdry` are generated.
The struct has four fields:
- `nel` --> number of elements associated with [`LineGpBasisFns`](@ref)
  `line_gp_fns`
- `uel_ids` --> mapping from `elem_id` to `unique_elem_id`
- `ufns` --> `num_unique_elems`×`ngp` matrix of unique 2-D basis functions, of
  type [`GpBasisFnsζα`](@ref)
- `bdry` --> [`Boundary`](@ref) of the parametric domain
"""
struct BdryGpBasisFns
  nel::Int64
  uel_ids::Vector{Int64}
  ufns::Matrix{GpBasisFnsζα}
  bdry::Boundary

  BdryGpBasisFns(
    line_gp_fns::LineGpBasisFns,
    perp_edge_fns::GpBasisFnsζ,
    bdry::Boundary
  ) = new(
    bdry_gp_basis_fns(line_gp_fns, perp_edge_fns, bdry)...
  );
end # BdryGpBasisFns


"""
    bdry_gp_basis_fns(line_gp_fns::LineGpBasisFns, perp_edge_fns::GpBasisFnsζ, bdry::Boundary)

External constructor of the `BdryGpBasisFns` `struct`.

Takes the tensor product of the line basis functions along the `bdry` boundary
with the basis functions along the perpendicular edge.

See also
[`BdryGpBasisFns`](@ref)
"""
function bdry_gp_basis_fns(
    line_gp_fns::LineGpBasisFns,
    perp_edge_fns::GpBasisFnsζ,
    bdry::Boundary
  )

  ufns    = line_gp_fns.ufns;
  num_uel = size(ufns, 1);
  num_gp  = size(ufns, 2);

  ugp_fns = Matrix{GpBasisFnsζα}(undef, num_uel, num_gp);
  for id = 1:num_uel, gp_id = 1:num_gp
    # account for which direction is ζ1 and which is ζ2
    ugp_fns[id, gp_id] = (bdry == BOTTOM || bdry == TOP) ?
      GpBasisFnsζα(ufns[id, gp_id], perp_edge_fns) :
      GpBasisFnsζα(perp_edge_fns, ufns[id, gp_id]);
  end

  return line_gp_fns.nel, line_gp_fns.uel_ids, ugp_fns, bdry;
end #bdry_gp_basis_fns


"""
    AreaGpBasisFns(line_gp_fns1::LineGpBasisFns, line_gp_fns2::LineGpBasisFns)

2-D B-spline basis functions and derivatives at all quadrature points in the
mesh area.

The two [`LineGpBasisFns`](@ref) structs passed in, `line_gp_fns1` and
`line_gp_fns2`, respectively contain 1-D derivatives in the ``ζ^1`` and ``ζ^2``
directions.
With the tensor product structure of B-splines over a surface
[piegl-tiller](@citep), it is straightforward to calculate the 2-D basis
functions and their derivatives.
As is the case for [`LineGpBasisFns`](@ref), a mapping from elements to unique
elements is generated; only unique basis function calculations are stored.
This struct has three fields:
- `nel` --> number of area elements
- `uel_ids` --> mapping from `elem_id` to `unique_elem_id`
- `ufns` --> `num_unique_elems`×`ngp` matrix of unique 2-D basis functions, of
  type [`GpBasisFnsζα`](@ref)
"""
struct AreaGpBasisFns
  nel::Int64
  uel_ids::Vector{Int64}
  ufns::Matrix{GpBasisFnsζα}

  AreaGpBasisFns(
    line_gp_fns1::LineGpBasisFns,
    line_gp_fns2::LineGpBasisFns
  ) = new(
    area_gp_basis_fns(line_gp_fns1, line_gp_fns2)...
  );
end # AreaGpBasisFns


"""
    area_gp_basis_fns(line_gp_fns1::LineGpBasisFns, line_gp_fns2::LineGpBasisFns)

External constructor of the `AreaGpBasisFns` struct.

Takes the tensor product of the boundary basis functions along the two
boundaries.

See also
[`AreaGpBasisFns`](@ref)
"""
function area_gp_basis_fns(
    line_gp_fns1::LineGpBasisFns,
    line_gp_fns2::LineGpBasisFns
  )

  num_el1  = line_gp_fns1.nel;      num_el2  = line_gp_fns2.nel;
  uel_ids1 = line_gp_fns1.uel_ids;  uel_ids2 = line_gp_fns2.uel_ids;
  ufns1    = line_gp_fns1.ufns;     ufns2    = line_gp_fns2.ufns;
  num_uel1 = size(ufns1, 1);        num_uel2 = size(ufns2, 1);
  num_gp1  = size(ufns1, 2);        num_gp2  = size(ufns2, 2);

  num_el  = num_el1 * num_el2;   # number of area elements
  num_gp  = num_gp1 * num_gp2;   # number of 2-D Gauss points per element
  num_uel = num_uel1 * num_uel2; # number of unique area elements


  uel_ids = zeros(Int64, num_el);
  for id1 = 1:num_el1, id2 = 1:num_el2
    uel_ids[id1 + (id2 - 1) * num_el1] =
      uel_ids1[id1] + (uel_ids2[id2] - 1) * num_uel1;
  end

  ugp_fns = Matrix{GpBasisFnsζα}(undef, num_uel, num_gp);
  for id1 = 1:num_uel1, id2 = 1:num_uel2
    for gp_id1 = 1:num_gp1, gp_id2 = 1:num_gp2
      ugp_fns[id1 + (id2 - 1) * num_uel1, gp_id1 + (gp_id2 - 1) * num_gp1] = 
        GpBasisFnsζα(ufns1[id1, gp_id1], ufns2[id2, gp_id2]);
    end
  end

  return num_el, uel_ids, ugp_fns;
end #area_gp_basis_fns


"""
    get_gp_basis_fn_weight(gpwξ::Float64, ζ::Float64, kv::KnotVector)::Float64

Returns the 1-D Gauss point weight scaled to the knot span length.

Integration by Gaussian quadrature is often carried out on the ξ line, with a
(dimensionless) length of 2. 
Here, integrals are evaluated on the knot span with length Δζ.
Thus, the weight in the ζ space is scaled by Δζ / 2.
"""
function get_gp_basis_fn_weight(
    gpwξ::Float64,
    ζ::Float64,
    kv::KnotVector,
  )::Float64

  index = get_knot_span_index(kv, ζ);
  Δζ = kv.ζs[index1+1] - kv.ζs[index1];

  return gpwξ * Δζ / 2;
end


"""
    get_gp_basis_fn_weight(gpwξ::Float64, ζ1::Float64, ζ2::Float64, kv1::KnotVector, kv2::KnotVector)::Float64

Returns the 2-D Gauss point weight scaled to the knot mesh area.

Integration by Gaussian quadrature is often carried out in the (ξ1, ξ2) square,
with a (dimensionless) area of 4. 
Here, integrals are evaluated in the knot mesh with area Δζ1 * Δζ2.
Thus, the weight in the (ζ1, ζ2) space is scaled by Δζ1 * Δζ2 / 4.
"""
function get_gp_basis_fn_weight(
    gpwξ::Float64,
    ζ1::Float64,
    ζ2::Float64, 
    kv1::KnotVector,
    kv2::KnotVector
  )::Float64
  index1 = get_knot_span_index(kv1, ζ1);
  index2 = get_knot_span_index(kv2, ζ2);

  Δζ1 = kv1.ζs[index1+1] - kv1.ζs[index1];
  Δζ2 = kv2.ζs[index2+1] - kv2.ζs[index2];

  return gpwξ * Δζ1 * Δζ2 / 4;
end


