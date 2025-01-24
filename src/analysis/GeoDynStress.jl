

export GeoDynStress;


"""
    GeoDynStress(xms, cps, dofs, N, ∂Nα, ∂∂Nαβ, p::Params)

Container for all geometric, dynamic, and stress-related quantities at a single
point ``(ζ^1, ζ^2)``.

Accessible fields:

|  code  | symbol |
|:------:|:------:|
`𝘅`      | ``\\bm{x}``
`𝗮_α`    | ``\\bm{a}_α``
`∂∂𝘅αβ`  | ``\\bm{x}_{, α β}``
`a_αβ`   | ``a_{α β}``
`a△αβ`   | ``a^{α β}``
`𝗮△α`    | ``\\bm{a}^α``
`JΩ`     | ``J_{\\Omega}``
`Γ△μ_αβ` | ``\\Gamma^\\mu_{\\alpha \\beta}``
`𝗻`      | ``\\bm{n}``
`b_αβ`   | ``b_{α β}``
`H`      | ``H``
`K`      | ``K``
`𝘃`      | ``\\bm{v}``
`∂𝘃α`    | ``\\bm{v}_{, α}``
`∂∂𝘃αβ`  | ``\\bm{v}_{, α β}``
`λ`      | ``λ``
`pm`     | ``p^{\\text{m}}``
`𝘃m`     | ``\\bm{v}^{\\text{m}}``
`∂𝘃mα`   | ``\\bm{v}^{\\text{m}}_{, α}``
`∂∂𝘃mαβ` | ``\\bm{v}^{\\text{m}}_{, α β}``
`𝞂`      | ``\\langle \\bm{\\sigma} \\rangle``
`𝞂m`     | ``\\langle \\bm{\\sigma}^{\\text{m}} \\rangle``
`𝗠`      | ``\\langle \\bm{M} \\rangle``
`𝗕a`     | ``[\\mathbf{B}^a]``
`𝗕b`     | ``[\\mathbf{B}^b]``
"""
struct GeoDynStress
  𝘅::SVector{XDIM,ComplexF64}
  𝗮_α::SMatrix{XDIM,ζDIM,ComplexF64}
  ∂∂𝘅αβ::SMatrix{XDIM,VOIGT,ComplexF64}
  a_αβ::SMatrix{ζDIM,ζDIM,ComplexF64}
  a△αβ::SMatrix{ζDIM,ζDIM,ComplexF64}
  𝗮△α::SMatrix{XDIM,ζDIM,ComplexF64}
  JΩ::ComplexF64
  Γ△μ_αβ::SMatrix{VOIGT,ζDIM,ComplexF64}
  𝗻::SVector{XDIM,ComplexF64}
  b_αβ::SMatrix{ζDIM,ζDIM,ComplexF64}
  H::ComplexF64
  K::ComplexF64
  𝘃::SVector{XDIM,ComplexF64}
  ∂𝘃α::SMatrix{XDIM,ζDIM,ComplexF64}
  ∂∂𝘃αβ::SMatrix{XDIM,VOIGT,ComplexF64}
  λ::ComplexF64
  pm::ComplexF64
  𝘃m::SVector{XDIM,ComplexF64}
  ∂𝘃mα::SMatrix{XDIM,ζDIM,ComplexF64}
  ∂∂𝘃mαβ::SMatrix{XDIM,VOIGT,ComplexF64}
  𝞂::SVector{VOIGT,ComplexF64}
  𝞂m::SVector{VOIGT,ComplexF64}
  𝗠::SVector{VOIGT,ComplexF64}
  𝗕a::Matrix{ComplexF64}
  𝗕b::Matrix{ComplexF64}

  GeoDynStress(
    xms::Matrix{ComplexF64},
    cps::Matrix{ComplexF64},
    dofs::Dict{Dof.Unknown, Int64},
    N::SVector{NEN,Float64},
    ∂Nα::SMatrix{NEN,ζDIM,Float64},
    ∂∂Nαβ::SMatrix{NEN,VOIGT,Float64},
    p::Params
  ) = new(
    geo_dyn_stress(xms, cps, dofs, N, ∂Nα, ∂∂Nαβ, p)...
  );
end # GeoDynStress


"""
    geo_dyn_stress(xms, cps, dofs, N, ∂Nα, ∂∂Nαβ, p::Params)

External constructor of the `GeoDynStress` `struct`.

See also
[`GeoDynStress`](@ref)
"""
function geo_dyn_stress(
    xms::Matrix{ComplexF64},
    cps::Matrix{ComplexF64},
    dofs::Dict{Dof.Unknown, Int64},
    N::SVector{NEN,Float64},
    ∂Nα::SMatrix{NEN,ζDIM,Float64},
    ∂∂Nαβ::SMatrix{NEN,VOIGT,Float64},
    p::Params
  )

  v_cps   = SMatrix{NEN,XDIM,ComplexF64}(
              get_dof_cps(cps, dofs, [Dof.vx, Dof.vy, Dof.vz]));
  vm_cps  = SMatrix{NEN,XDIM,ComplexF64}(
              get_dof_cps(cps, dofs, [Dof.vmx, Dof.vmy, Dof.vmz]));
  λ_cps   = SVector{NEN,ComplexF64}(get_dof_cps(cps, dofs, [Dof.λ]));
  p_cps   = SVector{NEN,ComplexF64}(get_dof_cps(cps, dofs, [Dof.pm]));


  # geometry

  𝘅       = SVector{XDIM,ComplexF64}(transpose(xms) * N);
  𝗮_α     = SMatrix{XDIM,ζDIM,ComplexF64}(transpose(xms) * ∂Nα);
  𝗮_1     = 𝗮_α[:,1];
  𝗮_2     = 𝗮_α[:,2];
  ∂∂𝘅αβ   = SMatrix{XDIM,VOIGT,ComplexF64}(transpose(xms) * ∂∂Nαβ);

  a_αβ    = transpose(𝗮_α) * 𝗮_α;
  a△αβ    = inv(a_αβ);
  𝗮△α     = 𝗮_α * a△αβ;
  JΩ      = sqrt(det(a_αβ));
  Γ△μ_αβ  = transpose(∂∂𝘅αβ) * 𝗮△α;

  𝗻       = cross(𝗮_1, 𝗮_2) / JΩ;
  b_αβ♭   = transpose(∂∂𝘅αβ) * 𝗻;
  b_αβ    = SMatrix{ζDIM,ζDIM,ComplexF64}(
              [b_αβ♭[1] b_αβ♭[3]; b_αβ♭[3] b_αβ♭[2]]);
  b△αβ    = a△αβ * b_αβ * a△αβ;

  H       = sum(a△αβ .* b_αβ) / 2;
  K       = det(b_αβ) / det(a_αβ);


  # dynamics

  𝘃       = SVector{XDIM,ComplexF64}(transpose(v_cps) * N);
  ∂𝘃α     = SMatrix{XDIM,ζDIM,ComplexF64}(transpose(v_cps) * ∂Nα);
  ∂∂𝘃αβ   = SMatrix{XDIM,VOIGT,ComplexF64}(transpose(v_cps) * ∂∂Nαβ);

  λ       = (transpose(λ_cps) * N)[1];
  pm      = (transpose(p_cps) * N)[1];

  𝘃m      = SVector{XDIM,ComplexF64}(transpose(vm_cps) * N);
  ∂𝘃mα    = SMatrix{XDIM,ζDIM,ComplexF64}(transpose(vm_cps) * ∂Nα);
  ∂∂𝘃mαβ  = SMatrix{XDIM,VOIGT,ComplexF64}(transpose(vm_cps) * ∂∂Nαβ);


  # stresses
  π△αβ    = transpose(𝗮△α) * (∂𝘃α * a△αβ);
  π△αβ    = (π△αβ + transpose(π△αβ)) * p.ζv;
  σ△αβ    = a△αβ * (p.kb*H^2 - p.kg*K + λ) - b△αβ * 2*p.kb*H + π△αβ;
  𝞂       = SVector{VOIGT,ComplexF64}([σ△αβ[1,1]; σ△αβ[2,2]; 2*σ△αβ[1,2]]);
  M△αβ    = a△αβ * H * (p.kb + 2*p.kg) - b△αβ * p.kg;
  𝗠       = SVector{VOIGT,ComplexF64}(
              [M△αβ[1,1]; M△αβ[2,2]; M△αβ[1,2] + M△αβ[2,1]]);

  σm△αβ   = transpose(𝗮△α) * (∂𝘃mα * a△αβ);
  σm△αβ   = (σm△αβ + transpose(σm△αβ)) * p.ζv;
  σm△αβ   = a△αβ * (p.kb*H^2 - p.kg*K) - b△αβ * 2*p.kb*H + σm△αβ;
  𝞂m      = SVector{VOIGT,ComplexF64}([σm△αβ[1,1]; σm△αβ[2,2]; 2*σm△αβ[1,2]]);


  # shape function matrices
  𝗕a      = zeros(ComplexF64, VOIGT, XDIM * size(∂Nα, 1));
  𝗕a[1,:] = reshape(𝗮_1 * transpose(∂Nα[:,1]), 1, :);
  𝗕a[2,:] = reshape(𝗮_2 * transpose(∂Nα[:,2]), 1, :);
  𝗕a[3,:] = reshape((
              𝗮_1 * transpose(∂Nα[:,2]) +
              𝗮_2 * transpose(∂Nα[:,1])
            )*0.5, 1, :);
  ∇∇Nαβ   = ∂∂Nαβ - ∂Nα * transpose(Γ△μ_αβ);
  𝗕b      = zeros(ComplexF64, VOIGT, XDIM * size(∂∂Nαβ, 1));
  for j=1:VOIGT
    𝗕b[j,:] = reshape(𝗻 * transpose(∇∇Nαβ[:,j]), 1, :);
  end
  
  return 𝘅, 𝗮_α, ∂∂𝘅αβ, a_αβ, a△αβ, 𝗮△α, JΩ, Γ△μ_αβ, 𝗻, b_αβ, H, K,
         𝘃, ∂𝘃α, ∂∂𝘃αβ, λ, pm, 𝘃m, ∂𝘃mα, ∂∂𝘃mαβ, 𝞂, 𝞂m, 𝗠, 𝗕a, 𝗕b
end # geo_dyn_stress


"""
    get_dof_cps(cps, dofs, dof_list)

Get control point values for `dof_list`, which may not be contained in `dofs`.

Allows one to evaluate various physical quantities (e.g. ``\\bm{v}``, ``λ``)
even if some of the associated control points are not degrees of freedom—and are
thus not solved for.
"""
function get_dof_cps(
    cps::Matrix{ComplexF64},
    dofs::Dict{Dof.Unknown, Int64},
    dof_list::Vector{Dof.Unknown}
  )::Matrix{ComplexF64}

  list_cps = zeros(ComplexF64, size(cps,1), length(dof_list));
  for (list_id, dof) in enumerate(dof_list)
    dof_id = get(dofs, dof, 0);
    if dof_id != 0
      list_cps[:,list_id] .= cps[:,dof_id];
    end
  end

  return list_cps;
end # get_dof_cps



