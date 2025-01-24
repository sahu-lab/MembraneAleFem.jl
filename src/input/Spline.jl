

export KnotVector


"""
    KnotVector(ζs::Vector{Float64}, nel::Int64, poly::Int64, curve::Curve)

Vector of knots, with which all 1-D B-spline functions are calculated.

The knot vector struct contains four fields:
- `ζs` --> list of all knots, including repeated ones
- `nel` --> number of 1-D elements, or knot spans
- `poly` --> polynomial order
- `curve` --> type of [Curve](@ref), which is either `CLAMPED` or `CLOSED`

All of these data are not required to generate a knot vector.
Two constructors are available—see [knot_vector](@ref)
"""
struct KnotVector
  ζs::Vector{Float64}
  nel::Int64
  poly::Int64
  curve::Curve
  
  KnotVector(ζs::Vector{Float64}, poly::Int64, curve::Curve) = new(
    knot_vector(ζs, poly, curve)...
  );

  KnotVector(nel::Int64, poly::Int64, curve::Curve) = new(
    knot_vector(nel, poly, curve)...
  );
end # KnotVector


"""
    knot_vector(ζs::Vector{Float64}, poly::Int64, curve::Curve)

External constructor of the [`KnotVector`](@ref) struct, where knots need not be
uniform.

The list of knots `ζs`, polynomial order `poly`, and [`Curve`](@ref) type
`curve` are passed in.
If they represent a valid knot vector, then the struct is created.
"""
function knot_vector(ζs::Vector{Float64}, poly::Int64, curve::Curve)

  if curve == CLAMPED
    @assert prod(ζs[1:end-1] .<= ζs[2:end]) "knots must be non-decreasing";
    @assert prod(ζs[1:poly+1] .== ζs[1]) "need poly+1 repeated knots at start";
    @assert prod(ζs[end-poly:end] .== ζs[end]) "need poly+1 repeated knots at end";
  elseif curve == CLOSED
    @assert ζs[2:end] - ζs[1:end-1] ≈ (ζs[2] - ζs[1]) * ones(length(ζs)-1);
  else
    @assert false "knot vector for $(curve) curve not implemented";
  end

  nel = length(ζs) -  2*poly - 1;

  return ζs, nel, poly, curve;
end # knot_vector


"""
    knot_vector(nel::Int64, poly::Int64, curve::Curve)

External constructor of the [`KnotVector`](@ref) struct with uniform knots.

The number of elements `nel`, polynomial order `poly`, and [`Curve`](@ref) type
`curve` are passed in.
If the `curve` is `CLOSED`, then no knots are repeated.
If the `curve` is `CLAMPED`, then the first and last `(poly+1)` knots
are repeated.
In both cases, the generated list of knots `ζs` ranges from 0 to 1.
"""
function knot_vector(nel::Int64, poly::Int64, curve::Curve)
  num_ζs = nel + 2 * poly + 1;
  ζs = zeros(Float64, num_ζs);

  if curve == CLAMPED
    # first and last (poly+1) knots are repeated
    for idx=poly+2:num_ζs-poly-1
      ζs[idx] = (idx-poly-1) / nel;
    end
    ζs[num_ζs-poly:num_ζs] .= 1.0;
  elseif curve == CLOSED
    # no repeated knots, with curve domain 0 ≤ ζ ≤ 1, for which the first poly
    # knots are < 0 and the last poly knots are > 1
    for idx=1:num_ζs
      ζs[idx] = (idx - poly - 1) / nel;
    end
  end

  # return by calling default constructor so assertions are checked as well
  return knot_vector(ζs, poly, curve);
end # knot_vector


# override == to compare each field
# from https://discourse.julialang.org/t/comparing-julia-structs/101554/2
Base.:(==)(kv1::KnotVector, kv2::KnotVector) =
  isequal(kv1.ζs, kv2.ζs) &&
  kv1.poly  == kv2.poly &&
  kv1.curve == kv2.curve;


# hash must also be altered, since == is used to determine if keys are equal
function Base.hash(kv::KnotVector, h_in::UInt)
  ht = hash(KnotVector, h_in);
  h1 = hash(kv.ζs, ht);
  h2 = hash(kv.poly, h1);
  h3 = hash(kv.curve, h2);
  return h3;
end


"""
    get_fine_ζs(nel::Int64, poly::Int64)::Vector{Float64}

Return knot list with knots concentrated in the center.
"""
function get_fine_ζs(
    nel::Int64,
    poly::Int64
  )::Vector{Float64}

  @assert nel >= 18 "fine mesh requires at least 18 1-D elements";

  num_ζs = nel + 2 * poly + 1;
  ζs = zeros(Float64, num_ζs);

  num_wide_l1 = 6;
  num_fine_l1 = nel - 2*num_wide_l1 - 1;

  ζ_wide_l1 = 1/3;
  ζ_fine_l1 = 1.0 - 2*ζ_wide_l1;

  
  # level 1, left portion
  for idx=poly+2:poly+num_wide_l1+1
    Δζ = ζ_wide_l1 / num_wide_l1;
    ζs[idx] = (idx-poly-1) * Δζ;
  end

  # level 1, right portion
  for idx=num_ζs-poly-num_wide_l1:num_ζs-poly-1
    Δζ = ζ_wide_l1 / num_wide_l1;
    ζs[idx] = 1.0 - ζ_wide_l1 + (idx - num_ζs + poly + num_wide_l1) * Δζ;
  end


  num_wide_l2 = floor(Int64, num_fine_l1 / 4);
  num_fine_l2 = num_fine_l1 - 2*num_wide_l2;

  ζ_wide_l2 = 1/9;
  ζ_fine_l2 = ζ_fine_l1 - 2*ζ_wide_l2;

  if 2*num_wide_l2 ≤ num_wide_l1
    # level 1, center portion
    Δζ = ζ_fine_l1 / (num_fine_l1 + 1);
    for idx = poly + num_wide_l1 + 2 : num_ζs - poly - num_wide_l1 - 1
      ζs[idx] = ζ_wide_l1 + (idx - poly - num_wide_l1 - 1) * Δζ;
    end
  else
    ## level 2, left portion
    Δζ = ζ_wide_l2 / num_wide_l2;
    for idx = num_wide_l1 + poly + 2 : num_wide_l1 + poly + 1 + num_wide_l2
      ζs[idx] = ζ_wide_l1 + (idx - num_wide_l1 - poly - 1) * Δζ;
    end
    ## level 2, right portion
    for idx = num_ζs - poly - num_wide_l1 - num_wide_l2 : num_ζs - poly - num_wide_l1 - 1
      ζs[idx] = 1.0 - ζ_wide_l1 - ζ_wide_l2 + (idx - num_ζs + poly + num_wide_l1 + num_wide_l2) * Δζ;
    end
    # level 2, center portion
    Δζ = ζ_fine_l2 / (num_fine_l2 + 1);
    for idx = num_wide_l1 + poly + 2 + num_wide_l2 : num_ζs - poly - num_wide_l1 - num_wide_l2 - 1
      ζs[idx] = ζ_wide_l1 + ζ_wide_l2 + (idx - num_wide_l1 - poly - 1 - num_wide_l2) * Δζ;
    end
  end

  ζs[num_ζs-poly:num_ζs] .= 1.0;
  return ζs;
end # get_fine_ζs


"""
    get_knot_span_index(kv::KnotVector, ζ::Float64)::Int64

Return knot span index of `ζ` in `CLAMPED` or `CLOSED` [`KnotVector`](@ref)
`kv`.

Algorithm A2.1 in Chap. 2, §5 of [piegl-tiller](@citet) ensures
`kv.ζs[idx]` ≤ `ζ` < `kv.ζs[idx+1]`,
where `idx` is the returned value.
The algorithm, valid only for `CLAMPED` knots, is here extended to `CLOSED`
knots as well.
"""
function get_knot_span_index(
    kv::KnotVector,
    ζ::Float64
  )::Int64

  # to adapt the CLAMPED algorithm to CLOSED knot vectors, we create a dummy
  # list of knots following the same format
  ζs = copy(kv.ζs);
  if kv.curve == CLOSED
    ζs[1:kv.poly] .= ζs[kv.poly + 1];
    ζs[end-kv.poly+1:end] .= ζs[end-kv.poly];
  end

  # for closed knot vectors, the first and last `poly` knots are excluded from
  # the domain; for clamped knot vectors these values are identical
  @assert ζ >= ζs[1]   "ζ smaller than first active knot";
  @assert ζ <= ζs[end] "ζ larger than last active knot";

  num_knots = length(ζs);
  # index of second-smallest knot
  m = 1;
  while ζs[m] == ζs[1]
    m = m+1;
  end

  # index of second-largest knot
  n = num_knots;
  while ζs[n] == ζs[end]
    n = n-1;
  end

  # check if ζ = ζs[end]
  if ζ == ζs[n+1]
    return n;
  end

  low  = m - 1;
  high = n + 1;

  # execute binary search
  mid  = floor(Int64, (low+high)/2);
  while ζ < ζs[mid] || ζ >= ζs[mid+1]
    if ζ < ζs[mid]
      high = mid;
    else
      low = mid;
    end
    mid  = floor(Int64, (low+high)/2);
  end
  return mid;
end # get_knot_span_index


"""
    get_bspline_vals(kv::KnotVector, ζ::Float64)::Vector{Float64}

Return `kv.poly`+1 nonzero B-spline functions at `ζ` for the given
[`KnotVector`](@ref) `kv`.

Algorithm A2.2 in Chap. 2, §5 of [piegl-tiller](@citet), valid for `CLAMPED`
knot vectors with (`kv.poly`+1) repeated knots at the start and end, is also
used for `CLOSED` knot vectors.
The local-to-global mapping returned by [`get_bspline_indices`](@ref) is
necessary to determine which global basis functions are nonzero at the provided
`ζ`.
"""
function get_bspline_vals(
    kv::KnotVector,
    ζ::Float64
  )::Vector{Float64}

  @assert ζ >= kv.ζs[1]   "ζ is less than smallest knot";
  @assert ζ <= kv.ζs[end] "ζ is greater than largest knot";

  poly  = kv.poly;
  ζs = kv.ζs;
  left  = zeros(Float64, poly + 1);
  right = zeros(Float64, poly + 1);

  span_idx = get_knot_span_index(kv, ζ);
  basis_functions = zeros(Float64, poly + 1);

  basis_functions[1] = 1.0;
  for j=1:poly
    left[j+1]  = ζ - ζs[span_idx + 1 - j];
    right[j+1] = ζs[span_idx + j] - ζ;
    saved    = 0.0;

    for r=1:j
      temp = basis_functions[r] / (right[r+1] + left[j+2-r]);
      basis_functions[r] = saved + right[r+1] * temp;
      saved = left[j+2-r] * temp;
    end

    basis_functions[j+1] = saved;
  end

  return basis_functions;
end # get_bspline_vals


"""
    get_bspline_ders(kv::KnotVector, ζ::Float64, num_ders::Int64)::Matrix{Float64}

Return B-spline functions and `num_ders` derivatives at `ζ` for the given
[`KnotVector`](@ref) `kv`.

Algorithm A2.3 in Chap. 2, §5 of [piegl-tiller](@citet), valid for `CLAMPED`
knot vectors with (`kv.poly`+1) repeated knots at the start and end, is also
used for `CLOSED` knot vectors.
Here `num_ders` is required to not be greater than `kv.poly`, and also
greater than or equal to zero.
The local-to-global mapping returned by [`get_bspline_indices`](@ref) is
necessary to determine which global basis functions are nonzero at the provided
`ζ`.

The returned matrix has `kv.poly`+1 rows and `num_ders`+1 columns.
The first column contains the functions themselves, and is thus identical to the
output of [`get_bspline_vals`](@ref).
The second column contains the first derivatives, and so on.
"""
function get_bspline_ders(
    kv::KnotVector,
    ζ::Float64,
    num_ders::Int64,
  )::Matrix{Float64}

  @assert ζ >= kv.ζs[1]   "ζ is less than smallest knot";
  @assert ζ <= kv.ζs[end] "ζ is greater than largest knot";

  poly  = kv.poly;
  ζs    = kv.ζs;
  left  = zeros(Float64, poly + 1);
  right = zeros(Float64, poly + 1);

  @assert num_ders >= 0    "cannot have fewer than zero derivatives";
  @assert num_ders <= poly "basis functions have only `poly` derivatives";

  span_idx = get_knot_span_index(kv, ζ);
  basis_fn_ders = zeros(Float64, poly + 1, num_ders + 1);

  # temporary matrix containing shape function values at polynomial orders
  ndu = zeros(Float64, poly + 1, poly + 1);
  
  # temporary matrix to store (alternating) most recently computed rows
  a = zeros(Float64, 2, poly + 1);


  ## compute B-spline values
  ndu[1,1] = 1.0;
  for j=1:poly
    left[j+1]  = ζ - ζs[span_idx + 1 - j];
    right[j+1] = ζs[span_idx + j] - ζ;
    saved    = 0.0;

    for r=1:j
      # lower triangle
      ndu[j+1,r] = right[r+1] + left[j+2-r];
      temp = ndu[r,j] / ndu[j+1,r];
      # upper triangle
      ndu[r,j+1] = saved + right[r+1] * temp;
      saved = left[j+2-r] * temp;
    end

    ndu[j+1,j+1] = saved;
  end


  ## load B-spline values
  for j=1:poly+1
    basis_fn_ders[j,1] = ndu[j,poly+1];
  end


  ## compute B-spline derivatives

  # loop over function index
  for r=1:poly+1

    # alternating rows
    s1 = 1; s2 = 2;
    a[1,1] = 1.0;

    # loop to compute k-th derivative
    for k=1:num_ders

      d = 0.0;
      rk = r - k;
      pk = poly - k;

      if r > k
        a[s2,1] = a[s1,1] / ndu[pk+2,rk];
        d = a[s2,1] * ndu[rk,pk+1];
      end

      j1 = rk >= 0 ? 1 : -rk + 1;
      j2 = r-2 <= pk ? k - 1 : poly - r + 1;

      for j=j1:j2
        a[s2,j+1] = (a[s1,j+1] - a[s1,j]) / ndu[pk+2,rk+j];
        d += a[s2,j+1] * ndu[rk+j,pk+1];
      end

      if r-1 <= pk
        a[s2,k+1] = -a[s1,k] / ndu[pk+2,r];
        d += a[s2,k+1] * ndu[r,pk+1];
      end

      basis_fn_ders[r,k+1] = d;

      # switch rows
      temp=s1; s1=s2; s2=temp;
    end
  end

  r = poly;
  for k=1:num_ders
    for j=1:poly+1
      basis_fn_ders[j,k+1] *= r;
    end
    r *= (poly-k);
  end

  return basis_fn_ders;
end # get_bspline_ders


"""
    get_bspline_indices(kv::KnotVector, ζ::Float64)::Vector{Int64}

Return global indices of nonzero B-splines at the value `ζ` in the
[`KnotVector`](@ref) `kv`.

By definition, there are only `kv.poly`+1 nonzero B-spline functions at any `ζ`.
The functions [`get_bspline_vals`](@ref) and [`get_bspline_ders`](@ref)
calculate these nonzero quantities, and the mapping returned here places them
within the global ordering.
The global indices are the same for all `ζ` within a single knot span, and thus
we begin by calling [`get_knot_span_index`](@ref).

Both `CLAMPED` knot vectors, with (`kv.poly`+1) repeats at the start and end, as
well as `CLOSED` knot vectors, are handled.
"""
function get_bspline_indices(
    kv::KnotVector,
    ζ::Float64
  )::Vector{Int64}

  ks_id = get_knot_span_index(kv, ζ);
  return get_bspline_indices(kv, ks_id);

end # get_bspline_indices


"""
    get_bspline_indices(kv::KnotVector, ks_id::Int64)::Vector{Int64}

Return global indices of nonzero B-splines at `ks_id```^{\\textrm{th}}`` knot
span in the [`KnotVector`](@ref) `kv`.

By definition, there are only `kv.poly`+1 nonzero B-spline functions over any
knot span.
The mapping returned here places these functions within the global ordering.
Both `CLAMPED` knot vectors, with (`kv.poly`+1) repeats at the start and end, as
well as `CLOSED` knot vectors, are handled.
Note that `ks_id` is *not* the element ID `eid`; rather, `ks_id = eid + kv.poly`
"""
function get_bspline_indices(
    kv::KnotVector,
    ks_id::Int64
  )::Vector{Int64}

  poly = kv.poly;
  ids = collect(1:poly+1) .+ (ks_id - poly - 1);

  if kv.curve == CLAMPED
    return ids;
  elseif kv.curve == CLOSED
    return (ids .- 1) .% kv.nel .+ 1;
  else
    @assert false "B-spline indices for $(kv.curve) curve not implemented";
    return Nothing;
  end
end # get_bspline_indices


"""
    get_1d_bspline_cps(kv::KnotVector, x::Function)::Vector{Float64}

Return the global control points ``\\{ x^{}_K \\}`` for the function ``x(ζ)``.

The global control points ``\\{ x^{}_K \\}`` are calculated for the provided
[`KnotVector`](@ref) `kv` such that
```math
x(ζ)
\\, \\approx \\, \\sum_{K = 1}^{\\texttt{nn}} N^{}_K (ζ) \\, x^{}_K
~,
```
where ``x(ζ)`` is the provided `x::Function` and here ``\\texttt{nn}`` is the
number of global 1-D nodes.
To determine the unknown ``\\{ x^{}_K \\}``, a set of ``\\texttt{nn}``
collocation points ``\\{ ζ^{}_J \\}`` is chosen via [`collocate_ζ`](@ref).
A matrix equation is then obtained as
```math
[x (ζ^{}_J)]
\\, = \\, [N^{}_K (ζ^{}_J)] \\, [x^{}_K]
~,
```
where the ``N^{}_K (ζ^{}_J)`` are determined with [`get_bspline_indices`](@ref)
and [`get_bspline_vals`](@ref).
The ``x^{}_K`` are then calculated and returned as a vector.
"""
function get_1d_bspline_cps(
    kv::KnotVector,
    x::Function
  )::Vector{Float64}

  ζlist     = collocate_ζ(kv);
  num_basis = length(ζlist);

  # basis function values evaluated at ζlist
  # TODO: use sparse matrix
  mat_basis = zeros(Float64, num_basis, num_basis);
  for j=1:num_basis # rows in mat_basis are at fixed ζ_j
    mat_basis[j, get_bspline_indices(kv, ζlist[j])] =
      get_bspline_vals(kv, ζlist[j]);
  end

  # matrix inverse for control points
  return mat_basis \ x.(ζlist);
end


"""
    get_2d_bspline_cps(kv1::KnotVector, kv2::KnotVector, x::Function)::Matrix{Float64}

Return the global control points ``\\{ x^{}_K \\}`` for the scalar-valued
function ``z(ζ^1, ζ^2)``.

The procedure of [`get_1d_bspline_cps`](@ref) is repeated, except in this case
collocation points are chosen across the 2-D parametric domain ``(ζ^1, ζ^2)``.
"""
function get_2d_bspline_cps(
    kv1::KnotVector,
    kv2::KnotVector,
    x::Function
  )::Matrix{Float64}

  ζ1list = collocate_ζ(kv1);
  ζ2list = collocate_ζ(kv2);
  num1   = length(ζ1list);
  num2   = length(ζ2list);
  num_basis = num1 * num2;

  # basis function values evaluated at ζ1list x ζ2list
  mat_basis = zeros(Float64, num_basis, num_basis);
  for k=1:num2, j=1:num1
    # relies on knowledge of the tensor product structure; not ideal
    indices = [i1 + (i2-1)*num1 for i2 in get_bspline_indices(kv2, ζ2list[k]),
               i1 in get_bspline_indices(kv1, ζ1list[j])];
    vals    = [v1*v2 for v2 in get_bspline_vals(kv2, ζ2list[k]),
               v1 in get_bspline_vals(kv1, ζ1list[j])];
    mat_basis[j + (k-1)*num1, indices[:,:]] = vals[:,:];
  end

  xvals = [x.(ζ1, ζ2) for ζ2 in ζ2list, ζ1 in ζ1list];

  # matrix inverse for control points
  return mat_basis \ reshape(transpose(xvals), :, 1);
end # get_2d_bspline_cps


"""
    collocate_ζ(kv::KnotVector)::Vector{Float64}

Take a [`KnotVector`](@ref) `kv` and return a list of points `ζ` used for
collocation.

The length of the returned list of ζ values is equal to the number of
global basis functions associated with the provided knot vector.
Both `CLAMPED` and `CLOSED` knot vectors are accounted for.
"""
function collocate_ζ(
    kv::KnotVector
  )::Vector{Float64}

  if kv.curve == CLOSED
    # closed curves have a simple collocation: the center of each knot span,
    # irrespective of the polynomial order
    return (kv.ζs[kv.poly+1:end-kv.poly-1] + kv.ζs[kv.poly+2:end-kv.poly]) / 2;

  elseif kv.curve == CLAMPED
    poly = kv.poly;
    ζs = unique(kv.ζs);

    @assert poly > 1 "interpolation for poly ≥ 2 only when clamped";
    @assert poly < 4 "interpolation for poly ≤ 3 only when clamped";
    @assert length(ζs) == length(kv.ζs) - 2*poly "no repeated interior knots";

    num_basis = length(kv.ζs) - poly - 1;
    ζlist     = zeros(Float64, num_basis);

    # for open knot vectors with poly+1 repeated first/last knots, basis
    # functions at the edge of the domain are interpolatory
    ζlist[1]   = ζs[1];
    ζlist[end] = ζs[end];

    if poly == 2
      # midpoint of each knot span
      for i=2:num_basis-1
        ζlist[i] = (ζs[i-1] + ζs[i]) / 2;
      end
    elseif poly == 3
      # midpoint of first and last knot span
      ζlist[2]     = (ζs[1] + ζs[2]) / 2;
      ζlist[end-1] = (ζs[end-1] + ζs[end]) / 2;
      # edge of each interior knot span
      for i=3:num_basis-2
        ζlist[i] = ζs[i-1];
      end
    end

    return ζlist;

  else
    @assert false "interpolation for $(kv.curve) curve not implemented";
    return Nothing;
  end
end # collocate_ζ


"""
    get_unique_1d_elements(kv::KnotVector)

Get all unique elements for a given knot vector `kv`.

This function returns several quantities, in the following order:
- `uel_num::Int64` -> number of unique elements
- `num_el::Int64` -> number of elements
- `uel_ids::Vector{Int64}` -> mapping from element id to unique element id
- `uel_list::Vector{Tuple{Float64, Float64}}` -> start and end `ζ` for each unique element
"""
function get_unique_1d_elements(
    kv::KnotVector
  )

  poly = kv.poly;
  num_knots = length(kv.ζs);
  prev_context = Tuple(zeros(Float64, 2*poly+1));

  num_el   = num_knots-2*poly-1;
  uel_ids  = zeros(Int64, num_el);
  uel_num  = 0;
  uel_list = Tuple{Float64, Float64}[];

  for k_idx=poly+1:num_knots-poly-1
    el_idx = k_idx - poly;
    # `context` contains the neighboring knot spans that affect basis functions
    context  = Tuple(kv.ζs[k_idx-poly+1 : k_idx+poly+1]
                   - kv.ζs[k_idx-poly   : k_idx+poly  ]);

    # '≈' accounts for floating point precision error
    if !prod(context .≈ prev_context)
      # new unique element
      uel_num += 1;
      push!(uel_list, (kv.ζs[k_idx], kv.ζs[k_idx+1]));
    end

    uel_ids[el_idx] = uel_num;
    prev_context = context;
  end

  return uel_num, num_el, uel_ids, uel_list
end # get_unique_1d_elements


"""
calc_spline_radius(kv::KnotVector, radius::Float64)::Float64

Return the control point radius corresponding to the physical `radius`.
"""
function calc_spline_radius(
    kv::KnotVector,
    radius::Float64
  )::Float64

  @assert kv.curve == CLOSED "can only calculate the radius for a closed curve";

  ζs, fns = get_global_bspline_ders(kv, 100, 0);
  xms = zeros(Float64, kv.nel, 2);
  for i = 1:kv.nel
    xms[i,1] = cos(pi*(1+2*i)/kv.nel);
    xms[i,2] = sin(pi*(1+2*i)/kv.nel);
  end

  xs = fns * xms;
  r_unit = sum(sqrt.(sum(xs.^2, dims=2))) / length(ζs);

  return radius / r_unit;
end # calc_spline_radius


"""
    get_global_bspline_ders(kv::KnotVector, num_ζ_intervals::Int64, d_order::Int64)

Calculate global B-splines and derivatives, including functions which are zero.

This function is used for plotting, as global basis functions are never used in
a numerical implementation.
Consequently, there are no unit tests for this function.

See also
[`KnotVector`](@ref),
[`get_bspline_ders`](@ref)
"""
function get_global_bspline_ders(
    kv::KnotVector,
    num_ζ_intervals::Int64,
    d_order::Int64
  )

  poly = kv.poly;

  if kv.curve == CLAMPED
    num_fns = length(kv.ζs) - 1 - poly;
    ζ_min = kv.ζs[1];
    ζ_max = kv.ζs[end];

  elseif kv.curve == CLOSED
    num_fns = kv.nel;
    ζ_min = kv.ζs[poly+1];
    ζ_max = kv.ζs[end-poly];

  else
    @assert false "bspline calculation for $(kv.curve) curve not implemented";
    return Nothing;
  end

  global_ders = zeros(Float64, num_ζ_intervals+1, num_fns);
  ζs = zeros(Float64, num_ζ_intervals+1);
  dζ = (ζ_max - ζ_min) / num_ζ_intervals;

  for ζ_idx = 1:num_ζ_intervals+1
    ζ = ζ_min + (ζ_idx - 1) * dζ;
    ζs[ζ_idx] = ζ;
    bspline_ders = get_bspline_ders(kv, ζ, poly);
    bspline_indices = get_bspline_indices(kv, ζ);
    global_ders[ζ_idx, bspline_indices] = bspline_ders[:, d_order+1];
  end

  return ζs, global_ders;
end # get_global_bspline_ders


