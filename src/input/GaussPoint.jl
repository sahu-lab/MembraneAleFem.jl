

export GaussPointsξ, GaussPointsζ


"""
    GaussPointsξ(ngp::Int64)

Weights ``\\hat{w}_k`` and points ``ξ_k`` on the interval ``[-1, +1]`` with
which a 1-D integral is approximated.

The struct has three fields:
- `ngp` --> number of Gaussian quadrature points
- `ξs`  --> set of points ``ξ_k``
- `ws`  --> set of weights ``\\hat{w}_k``
Calculations from [Wikipedia: Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Legendre_quadrature).
"""
struct GaussPointsξ
  ngp::Int64
  ξs::Vector{Float64}
  ws::Vector{Float64}

  GaussPointsξ(ngp::Int64) = new(
    ngp,
    get_ξs(ngp),
    get_ξ_weights(ngp),
  );

end # GaussPointsξ


# override == to compare each field
# from https://discourse.julialang.org/t/comparing-julia-structs/101554/2
Base.:(==)(gpξ1::GaussPointsξ, gpξ2::GaussPointsξ) =
  gpξ1.ngp == gpξ2.ngp &&
  isequal(gpξ1.ξs, gpξ2.ξs) &&
  isequal(gpξ1.ws, gpξ2.ws);


# hash must also be altered, since == is used to determine if keys are equal
function Base.hash(gpξ::GaussPointsξ, h_in::UInt)
  ht = hash(GaussPointsξ, h_in);
  h1 = hash(gpξ.ngp, ht);
  h2 = hash(gpξ.ξs,  h1);
  h3 = hash(gpξ.ws,  h2);
  return h3;
end


"""
    GaussPointsζ(gpsξ::GaussPointsξ, ζlo::Float64, ζhi::Float64)

Weights ``w_k`` and points ``ζ_k`` with which a 1-D integral is approximated
over `ζlo` ≤ `ζ` ≤ `ζhi`.

Calculated via the transform from `ξ` to `ζ` via the relations

`ζ```_k`` = `ξ```_k`` (`ζhi` - `ζlo`) / 2 + (`ζhi` + `ζlo`) / 2

``w_k`` = ``\\hat{w}_k`` (`ζhi` - `ζlo`) / 2

The struct has three fields:
- `ngp` --> number of Gaussian quadrature points
- `ζs`  --> set of points ``ζ_k``
- `ws`  --> set of weights ``w_k``
"""
struct GaussPointsζ
  ngp::Int64
  ζs::Vector{Float64}
  ws::Vector{Float64}

  GaussPointsζ(gpsξ::GaussPointsξ, ζlo::Float64, ζhi::Float64) = new(
    gpsξ.ngp,
    gpsξ.ξs * (ζhi - ζlo) / 2 .+ (ζhi + ζlo) / 2,
    gpsξ.ws * (ζhi - ζlo) / 2
  );

end # GaussPointsζ


# override == to compare each field
# from https://discourse.julialang.org/t/comparing-julia-structs/101554/2
Base.:(==)(gpζ1::GaussPointsζ, gpζ2::GaussPointsζ) =
  gpζ1.ngp == gpζ2.ngp &&
  isequal(gpζ1.ζs, gpζ2.ζs) &&
  isequal(gpζ1.ws, gpζ2.ws);


# hash must also be altered, since == is used to determine if keys are equal
function Base.hash(gpζ::GaussPointsζ, h_in::UInt)
  ht = hash(GaussPointsζ, h_in);
  h1 = hash(gpζ.ngp, ht);
  h2 = hash(gpζ.ζs,  h1);
  h3 = hash(gpζ.ws,  h2);
  return h3;
end


"""
    get_ξs(ngp::Int64)

Returns a vector of `ngp` Gaussian quadrature points on the domain -1 ≤ `ξ` ≤ 1.

Calculations from [Wikipedia: Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Legendre_quadrature).
"""
function get_ξs(
    ngp::Int64
  )::Vector{Float64}

  @assert ngp <= 4 "have at most 4 1-D Gauss points";

  if ngp == 1
    # 1 1-D point, for generation of element grid only
    return [ -1.0 ];

  elseif ngp == 3
    # 3 1-D points
    return [ -sqrt(3/5), 0.0, sqrt(3/5) ];

  elseif ngp == 4
    # 4 1-D points
    return [
            -sqrt(3/7 + sqrt(6/5)*2/7),
            -sqrt(3/7 - sqrt(6/5)*2/7),
             sqrt(3/7 - sqrt(6/5)*2/7),
             sqrt(3/7 + sqrt(6/5)*2/7)
           ];
  else
    @assert false "ξs for $(ngp) Gauss points not implemented";
  end

  return Nothing;
end # get_ξs


"""
    get_ξ_weights(ngp::Int64)

Returns a vector of weights on the line -1 ≤ `ξ` ≤ 1.

Calculations from [Wikipedia: Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Legendre_quadrature).
"""
function get_ξ_weights(
    ngp::Int64
  )::Vector{Float64}

  @assert ngp <= 4 "have at most 4 1-D Gauss points";

  if ngp == 1
    # 1 1-D point, for generation of element grid only
    # the weight is set to 0, as this should not be used for integration
    return [ 0.0 ];

  elseif ngp == 3
    # 3 1-D points
    return [ 5/9, 8/9, 5/9];

  elseif ngp == 4
    # 4 1-D points
    return [
            (18 - sqrt(30))/36,
            (18 + sqrt(30))/36,
            (18 + sqrt(30))/36,
            (18 - sqrt(30))/36
           ];
  else
    @assert false "Weights for $(ngp) Gauss points not implemented";
  end
  
  return Nothing;
end # get_ξ_weights


