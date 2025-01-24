

"""
    time_step!(mesh, xms, cps, Δt, p::Params; args...)

Given a known state of the membrane, determine the state a time `Δt` later.

The Newton--Raphson method is employed to advance forward in time, which uses
the result of [`calc_r_K`](@ref MembraneAleFem.calc_r_K).
"""
function time_step!(
    mesh::Mesh,
    xms::Matrix{Float64},
    cps::Matrix{Float64},
    time::Float64,
    Δt::Float64,
    p::Params;
    args...
  )

  if p.output
    io = open(args[:out_path] * "/" * args[:out_file], "a");
  end

  εs = Float64[];
  iter = 1;

  # Newton--Raphson loop
  while iter < 15
    if p.output
      print(io, "--> iteration ", iter, ": ");
    end

    # calculate global residual vector and tangent diffusion matrix
    r_gl, K_gl = calc_r_K(mesh, xms, cps, time, Δt, p; args...);

    # solve for change in unknowns
    Δu = -K_gl \ r_gl;

    # change in unknown vector --> change in control point matrix
    Δcps = zeros(size(cps));
    Δcps[mesh.ID_inv[:]] .= Δu;

    # update control points and mesh position
    cps .+= Δcps;
    update_xms!(p.motion, xms, Δcps, Δt, mesh.dofs);

    push!(εs, norm(Δu)/mesh.nmdf);
    if p.output
      println(io, "ε = ", εs[end]);
    end

    iter += 1;
    if εs[end] < p.εnr break; end
  end # Newton--Raphson iteration

  if p.output
    close(io);
  end

  @assert εs[end] < p.εnr "did not reach Newton--Raphson tolerance";
  return Nothing;
end # time_step!


"""
    calc_r_K(mesh, xms, cps, time, Δt, p::Params; args...)

Calculate the residual vector r and diffusion matrix K.

The global tangent diffusion matrix K is calculated via numerical
differentiation, in which the residual vector and its argument are extended into
the complex plane [lyness-siam-1967](@citet), [lyness-mc-1968](@citet).
"""
function calc_r_K(
    mesh::Mesh,
    xms::Matrix{Float64},
    cps::Matrix{Float64},
    time::Float64,
    Δt::Float64,
    p::Params;
    args...
  )

  m_motion_order = get_m_motion_order(p.motion, mesh.dofs);

  # from https://julialang.org/blog/2023/07/PSA-dont-use-threadid
  chunk_size  = ceil(Int64, max(1, mesh.numel/nthreads()));
  data_chunks = partition(1:mesh.numel, chunk_size);

  tasks = map(data_chunks) do chunk
    @spawn begin

      r_th =   zeros(Float64, mesh.nmdf);
      K_th = spzeros(Float64, mesh.nmdf, mesh.nmdf);

      # loop over area elements
      for el_id ∈ chunk

        xms_el = complex(xms[mesh.IX[:,el_id], :]);
        cps_el = complex(cps[mesh.IX[:,el_id], :]);

        r_el = real(calc_elem_residual(mesh, el_id, xms_el, cps_el, p));
        K_el = zeros(Float64, p.nen * mesh.ndf, p.nen * mesh.ndf);

        # element active unknown IDs, with LM[:] giving the global unknown IDs
        el_unk_ids = findall(!=(0), mesh.LM[:, el_id]);

        # calculate element tangent diffusion matrix
        for (el_node_id, gl_node_id) ∈ enumerate(mesh.IX[:,el_id])
          for dof_id ∈ findall(!=(0), mesh.ID[:,gl_node_id])

            Δcps = copy(cps_el); Δcps[el_node_id, dof_id] += p.εk * im;
            K_el[:, dof_id+mesh.ndf*(el_node_id-1)] +=
              imag(calc_elem_residual(mesh, el_id, xms_el, Δcps, p)) / p.εk;

            # mesh velocity perturbations also perturb positions
            if dof_id ∈ m_motion_order
              Δxms = copy(xms_el);
              Δxms[el_node_id, findfirst(==(dof_id), m_motion_order)] += p.εk * im;
              K_el[:, dof_id+mesh.ndf*(el_node_id-1)] += (
                imag(calc_elem_residual(mesh, el_id, Δxms, cps_el, p)) / p.εk) * Δt;
            end

          end
        end

        # write to global residual and tangent
        for el_unk_rid ∈ el_unk_ids
          gl_unk_rid = mesh.LM[el_unk_rid, el_id];
          r_th[gl_unk_rid] += r_el[el_unk_rid];
          for el_unk_cid ∈ el_unk_ids
            gl_unk_cid = mesh.LM[el_unk_cid, el_id];
            K_th[gl_unk_rid, gl_unk_cid] += K_el[el_unk_rid, el_unk_cid];
          end
        end

      end # area element loop

      return r_th, K_th;
    end
  end

  th_rKs = fetch.(tasks);

  r_gl = sum(th_rKs[i][1] for i = 1:nthreads());
  K_gl = sum(th_rKs[i][2] for i = 1:nthreads());


  # Neumann boundary conditions
  for (bdry, ntype, nval) ∈ mesh.inh_neu_bcs

    # loop over boundary elements
    for el_id ∈ mesh.bdry_elems[bdry]

      xms_el = complex(xms[mesh.IX[:,el_id], :]);
      cps_el = complex(cps[mesh.IX[:,el_id], :]);

      r_el = real(calc_bdry_element_residual(mesh, bdry, ntype, nval, el_id,
                                             xms_el, cps_el, time, p; args...));
      K_el = zeros(Float64, p.nen * mesh.ndf, p.nen * mesh.ndf);

      # element active unknown IDs, with LM[:] giving the global unknown IDs
      el_unk_ids = findall(!=(0), mesh.LM[:, el_id]);

      # calculate element tangent diffusion matrix
      for (el_node_id, gl_node_id) ∈ enumerate(mesh.IX[:,el_id])
        for dof_id ∈ findall(!=(0), mesh.ID[:,gl_node_id])

          Δcps = copy(cps_el); Δcps[el_node_id, dof_id] += p.εk * im;
          K_el[:, dof_id+mesh.ndf*(el_node_id-1)] += imag(
            calc_bdry_element_residual(mesh, bdry, ntype, nval, el_id,
                                       xms_el, Δcps, time, p; args...)) / p.εk;

          # mesh velocity perturbations also perturb positions
          if dof_id ∈ m_motion_order
            Δxms = copy(xms_el);
            Δxms[el_node_id, findfirst(==(dof_id), m_motion_order)] += p.εk*im;
            K_el[:, dof_id+mesh.ndf*(el_node_id-1)] += (imag(
              calc_bdry_element_residual(mesh, bdry, ntype, nval, el_id, Δxms,
                                         cps_el, time, p; args...)) / p.εk)*Δt;
          end
        end
      end

      # write to global residual and tangent
      for el_unk_rid ∈ el_unk_ids
        gl_unk_rid = mesh.LM[el_unk_rid, el_id];
        r_gl[gl_unk_rid] += r_el[el_unk_rid];
        for el_unk_cid ∈ el_unk_ids
          gl_unk_cid = mesh.LM[el_unk_cid, el_id];
          K_gl[gl_unk_rid, gl_unk_cid] += K_el[el_unk_rid, el_unk_cid];
        end
      end

    end # boundary element loop
  end # loop over edges

  return r_gl, K_gl;
end # calc_r_K


"""
    calc_elem_residual(mesh, el_id, xms_el, cps_el, p::Params)

Calculate the area element residual vector by looping over Gauss points.
"""
function calc_elem_residual(
    mesh::Mesh,
    el_id::Int64,
    xms_el::Matrix{ComplexF64},
    cps_el::Matrix{ComplexF64},
    p::Params
  )::Vector{ComplexF64}

  λ_order = get_λ_order(mesh.dofs);
  v_order = get_v_order(mesh.dofs);
  m_order = get_m_order(mesh.dofs);
  p_order = get_p_order(mesh.dofs);

  rv_el, rm_el, rλ_el, rp_el = calc_elem_dof_residuals(mesh, el_id,
                                                       xms_el, cps_el, p);

  # construct local residual vector, according to mesh.dofs
  r_el = zeros(ComplexF64, p.nen * mesh.ndf);
  for el_node=1:p.nen
    for vj ∈ findall(!=(0), v_order)
      r_el[v_order[vj] + mesh.ndf * (el_node-1)] =
      rv_el[vj +  length(v_order) * (el_node-1)];
    end
    for mj ∈ findall(!=(0), m_order)
      r_el[m_order[mj] + mesh.ndf * (el_node-1)] =
      rm_el[mj +  length(m_order) * (el_node-1)];
    end
    r_el[λ_order + mesh.ndf * (el_node-1)] = rλ_el[el_node];
    if p_order != 0
      r_el[p_order + mesh.ndf * (el_node-1)] = rp_el[el_node];
    end
  end

  return r_el;
end # calc_elem_residual


"""
    calc_elem_dof_residuals(mesh, el_id, xms_el, cps_el, p::Params)

Calculate the area element residual vectors for each degree of freedom by looping over Gauss points.

The discretized form of the residual vectors can be read directly from this
function.
"""
function calc_elem_dof_residuals(
    mesh::Mesh,
    el_id::Int64,
    xms_el::Matrix{ComplexF64},
    cps_el::Matrix{ComplexF64},
    p::Params
  )

  λ_order = get_λ_order(mesh.dofs);
  v_order = get_v_order(mesh.dofs);
  m_order = get_m_order(mesh.dofs);
  p_order = get_p_order(mesh.dofs);

  # number of element degrees of freedom for unknowns
  nevdf = p.nen * length(v_order);
  nemdf = p.nen * length(m_order);
  neλdf = p.nen;
  nepdf = p_order == 0 ? 0 : p.nen;

  # elemental residual vectors
  rv_el = zeros(ComplexF64, nevdf);
  rm_el = zeros(ComplexF64, nemdf);
  rλ_el = zeros(ComplexF64, neλdf);
  rp_el = zeros(ComplexF64, nepdf);

  # Dohrmann--Bochev matrices
  GDB_el = @SMatrix zeros(ComplexF64, NEDBDF, NEN);
  HDB_el = @SMatrix zeros(ComplexF64, NEDBDF, NEDBDF);
  ξs_1d  = GaussPointsξ(GP1D).ξs;

  # loop over Gauss points
  for gp_id = 1:GP1D^ζDIM

    N     = get_N(el_id, gp_id, mesh);
    ∂Nα   = get_∂Nα(el_id, gp_id, mesh);
    ∂∂Nαβ = get_∂∂Nαβ(el_id, gp_id, mesh);
    gpw   = get_gpw(el_id, gp_id, mesh);

    gds = GeoDynStress(xms_el, cps_el, mesh.dofs, N, ∂Nα, ∂∂Nαβ, p);

    rv_el += ( transpose(gds.𝗕a) * gds.𝞂 +
               transpose(gds.𝗕b) * gds.𝗠 ) * gds.JΩ * gpw;
    if p.pn != 0.0
      rv_el -= reshape((gds.𝗻) * transpose(N), :, 1) * p.pn * gds.JΩ * gpw;
    end
    rλ_el += N * gpw * ( gds.JΩ * tr(transpose(gds.𝗮△α) * gds.∂𝘃α)
                         - p.αdb * gds.λ / p.ζv );
    if p.motion == EUL
      # the "dot" (⋅) operator involves a complex conjugate, so cannot be used
      𝘃_𝗻    = gds.𝗻 * sum(gds.𝗻 .* gds.𝘃);
      rm_el += reshape((gds.𝘃m - 𝘃_𝗻)*transpose(N), :, 1) * p.αm * gds.JΩ * gpw;
    elseif p.motion == ALEV || p.motion == ALEVB
      rm_el += ( transpose(gds.𝗕a) * gds.𝞂m ) * gds.JΩ * gpw;
      if p.motion == ALEVB
        rm_el += ( transpose(gds.𝗕b) * gds.𝗠  ) * gds.JΩ * gpw;
      end
      rm_el -= reshape((gds.𝗻) * transpose(N), :, 1) * gds.pm * gds.JΩ * gpw;
      # the "dot" (⋅) operator involves a complex conjugate, so cannot be used
      rp_el -= N * gpw * gds.JΩ * (sum(gds.𝗻 .* (gds.𝘃m - gds.𝘃)));
      rp_el -= N * gpw * p.αdb * gds.pm / p.ζv;
    end

    NDB = @SVector [ξs_1d[(gp_id - 1) % GP1D + 1];
                    ξs_1d[floor(Int64, (gp_id - 1) / GP1D) + 1];
                    1.0];

    GDB_el += NDB * transpose(N)   * gpw;
    HDB_el += NDB * transpose(NDB) * gpw;
  end # Gauss point loop

  tmpDB  = transpose(GDB_el) * inv(HDB_el) * GDB_el;
  rλ_el += tmpDB * cps_el[:, λ_order] * p.αdb / p.ζv;
  if p.motion == ALEV || p.motion == ALEVB
    rp_el += tmpDB * cps_el[:, p_order] * p.αdb / p.ζv;
  end

  return rv_el, rm_el, rλ_el, rp_el;
end # calc_elem_dof_residuals


"""
    calc_bdry_element_residual(mesh, bdry, ntype, nval, ..., p::Params; args...)

Calculate the boundary element residual vector by looping over Gauss points.
"""
function calc_bdry_element_residual(
    mesh::Mesh,
    bdry::Boundary,
    ntype::Neumann,
    nval::Float64,
    el_id::Int64,
    xms_el::Matrix{ComplexF64},
    cps_el::Matrix{ComplexF64},
    time::Float64,
    p::Params;
    args...
  )::Vector{ComplexF64}

  v_order = get_v_order(mesh.dofs);
  m_order = get_m_order(mesh.dofs);

  # number of element degrees of freedom for velocity, mesh velocity
  nevdf = NEN * length(v_order);
  nemdf = NEN * length(m_order);

  # elemental residual vectors
  rv_el = zeros(ComplexF64, nevdf);
  rm_el = zeros(ComplexF64, nemdf);

  # loop over Gauss points
  for gp_id = 1:GP1D

    N     = get_N(bdry, el_id, gp_id, mesh);
    ∂Nα   = get_∂Nα(bdry, el_id, gp_id, mesh);
    ∂∂Nαβ = get_∂∂Nαβ(bdry, el_id, gp_id, mesh);
    gpw   = get_gpw(bdry, el_id, gp_id, mesh);

    gds   = GeoDynStress(xms_el, cps_el, mesh.dofs, N, ∂Nα, ∂∂Nαβ, p);
    (τ,ν) = calc_τ_ν(bdry, gds.𝗮_α, gds.𝗻);
    JΓ    = 1 / sqrt(sum( (transpose(gds.𝗮△α) * τ).^2 ));

    if ntype == STRETCH || ntype == SHEAR
      𝗳 = ntype == STRETCH ? nval * ν : nval * τ;
      rv_el -= reshape(𝗳 * transpose(N), :, 1) * JΓ * gpw;
    elseif ntype == MOMENT && p.scenario == F_BEND
      ν△α  = transpose(gds.𝗮△α) * ν;
      Mval = nval * min(time/args[:bend_tm], 1.0);
      rv_el -= reshape(gds.𝗻 * transpose(∂Nα * ν△α), :, 1) * Mval * JΓ * gpw;
    else
      @assert false "Neumann boundary condition not implemented";
    end
  end # Gauss point loop

  # construct local residual vector, according to mesh.dofs
  r_el = zeros(ComplexF64, NEN * mesh.ndf);
  for el_node=1:NEN
    for vj ∈ findall(!=(0), v_order)
      r_el[v_order[vj] + mesh.ndf * (el_node-1)] =
      rv_el[vj +  length(v_order) * (el_node-1)];
    end
    for mj ∈ findall(!=(0), m_order)
      r_el[m_order[mj] + mesh.ndf * (el_node-1)] =
      rm_el[mj +  length(m_order) * (el_node-1)];
    end
  end

  return r_el;
end # calc_bdry_element_residual


"""
    update_xms!(motion, xms, cps, Δt, dofs)

Update membrane positions according to the mesh velocity, using forward Euler.
"""
function update_xms!(
    motion::Motion,
    xms::Matrix{Float64},
    cps::Matrix{Float64},
    Δt::Float64,
    dofs::Dict{Dof.Unknown, Int64}
  )

  m_motion_order = get_m_motion_order(motion, dofs);

  for mj ∈ findall(!=(0), m_motion_order)
    xms[:, mj] .+= Δt * cps[:, m_motion_order[mj]];
  end

  return;
end # update_xms!


"""
    calc_τ_ν(bdry::Boundary, 𝗮_α::Matrix{ComplexF64}, 𝗻::Vector{ComplexF64})

Determine the in-plane unit vectors `τ` and `ν` on the boundary.
"""
function calc_τ_ν(
    bdry::Boundary,
    𝗮_α::SMatrix{XDIM,ζDIM,ComplexF64},
    𝗻::SVector{XDIM,ComplexF64}
  )::Tuple{SVector{XDIM,ComplexF64}, SVector{XDIM,ComplexF64}}

  #τ = zeros(ComplexF64, XDIM);
  if bdry == BOTTOM
    τ = 𝗮_α[:,1];
  elseif bdry == RIGHT
    τ = 𝗮_α[:,2];
  elseif bdry == TOP
    τ = -𝗮_α[:,1];
  elseif bdry == LEFT
    τ = -𝗮_α[:,2];
  end

  τ /= sqrt(dot(τ, τ));
  ν = cross(τ, 𝗻);

  return (τ, ν);
end # calc_τ_ν


