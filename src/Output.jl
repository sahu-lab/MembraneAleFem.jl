

export generate_state, generate_output, analyze_output


"""
    generate_state(path::String, t_id::Int64)

Generate the membrane state from the given path and time id.
"""
function generate_state(
    path::String,
    t_id::Int64
  )

  p = deserialize("$(path)/params.dat");
  args = deserialize("$(path)/args.dat");

  mesh = MembraneAleFem.Mesh(p; args...);
  xms  = readdlm("$(path)/t$(t_id)-xms.txt");
  cps  = readdlm("$(path)/t$(t_id)-cps.txt");

  return mesh, xms, cps;
end # generate_state


"""
    generate_output(mesh::Mesh, xms::Matrix{Float64}, cps::Matrix{Float64})

Generate visualization output for the given mesh, positions, and unknowns.
"""
function generate_output(
    mesh::Mesh,
    xms::Matrix{Float64},
    cps::Matrix{Float64}
  )

  numζ1s = mesh.num1el * GP1D + 2;
  numζ2s = mesh.num2el * GP1D + 2;

  xout = zeros(Float64, numζ1s, numζ2s, XDIM);
  uout = zeros(Float64, numζ1s, numζ2s, mesh.ndf);

  # area Gauss points
  for el2_id=1:mesh.num2el, el1_id=1:mesh.num1el
    el_id = el1_id + mesh.num1el * (el2_id - 1);
    for gp2_id = 1:GP1D, gp1_id = 1:GP1D
      gp_id = gp1_id + GP1D * (gp2_id - 1);
      x = xms[mesh.IX[:,el_id], :]' * get_N(el_id, gp_id, mesh);
      u = cps[mesh.IX[:,el_id], :]' * get_N(el_id, gp_id, mesh);
      out1_id = (el1_id-1)*GP1D + gp1_id + 1;
      out2_id = (el2_id-1)*GP1D + gp2_id + 1;
      xout[out1_id, out2_id, :] = x;
      uout[out1_id, out2_id, :] = u;
    end
  end

  # boundary Gauss points
  for bdry in instances(Boundary)
    out_id = 2;
    for el_id in mesh.bdry_elems[bdry], gp_id = 1:GP1D
      x = xms[mesh.IX[:,el_id], :]' * get_N(bdry, el_id, gp_id, mesh);
      u = cps[mesh.IX[:,el_id], :]' * get_N(bdry, el_id, gp_id, mesh);

      # map from boundary to output arrays
      out1_id = 0; out2_id = 0;

      if bdry == BOTTOM
        out1_id = out_id;
        out2_id = 1;
      elseif bdry == RIGHT
        out1_id = numζ1s;
        out2_id = out_id;
      elseif bdry == TOP
        out1_id = out_id;
        out2_id = numζ2s;
      elseif bdry == LEFT
        out1_id = 1;
        out2_id = out_id;
      end

      # write to output arrays
      xout[out1_id, out2_id, :] = x;
      uout[out1_id, out2_id, :] = u;

      out_id += 1;
    end
  end

  # corner Gauss points
  for crnr in instances(Corner)
    el_id = mesh.crnr_elems[crnr];
    x = xms[mesh.IX[:,el_id], :]' * mesh.crnr_gp_fns[crnr].N;
    u = cps[mesh.IX[:,el_id], :]' * mesh.crnr_gp_fns[crnr].N;

    # map from corner to output arrays
    out1_id = 0; out2_id = 0;
    if crnr == BOTTOM_LEFT
      out1_id = 1;
      out2_id = 1;
    elseif crnr == BOTTOM_RIGHT
      out1_id = numζ1s;
      out2_id = 1;
    elseif crnr == TOP_LEFT
      out1_id = 1;
      out2_id = numζ2s;
    elseif crnr == TOP_RIGHT
      out1_id = numζ1s;
      out2_id = numζ2s;
    end

    # write to output arrays
    xout[out1_id, out2_id, :] = x;
    uout[out1_id, out2_id, :] = u;
  end

  return xout, uout;
end # generate_output


"""
    analyze_output(path::String, t_id::Int64)

Return user-specified information, with which to analyze the chosen data.

Currently, we do not have separate functions for different scenarios, and so it
is up to the user to decide how they would like to analyze their data.
"""
function analyze_output(
    path::String,
    t_id::Int64
  )

  p = deserialize("$(path)/params.dat");
  mesh, xms, cps = generate_state(path, t_id);

  iod = open(path * "/data.txt", "w");
  println(iod, "x\ty\tz\tr\tλ\tH\tK");

  # error analysis
  Er = 0.0; Eλ = 0.0; EH = 0.0;
  λ_ex = 0.25; H_ex = 0.5;

  for el2_id=1:mesh.num2el, el1_id=1:mesh.num1el
    el_id = el1_id + mesh.num1el * (el2_id - 1);
    xms_el = complex(xms[mesh.IX[:,el_id], :]);
    cps_el = complex(cps[mesh.IX[:,el_id], :]);

    for gp2_id = 1:GP1D, gp1_id = 1:GP1D
      gp_id = gp1_id + GP1D * (gp2_id - 1);

      N     = get_N(el_id, gp_id, mesh);
      ∂Nα   = get_∂Nα(el_id, gp_id, mesh);
      ∂∂Nαβ = get_∂∂Nαβ(el_id, gp_id, mesh);
      gpw   = get_gpw(el_id, gp_id, mesh);

      gds = GeoDynStress(xms_el, cps_el, mesh.dofs, N, ∂Nα, ∂∂Nαβ, p);

      𝘅 = real(gds.𝘅);  𝘃 = real(gds.𝘃);  𝘃m = real(gds.𝘃m);
      λ = real(gds.λ);  H = real(gds.H);   K = real(gds.K);
      𝗻 = real(gds.𝗻);

      Er += gpw * ( √((𝘅[1]-sin(1))^2+(𝘅[3]-cos(1))^2) - 1.0)^2
      Eλ += gpw * (λ - λ_ex)^2;
      EH += gpw * (H - H_ex)^2;

      println(iod, 𝘅[1], "\t", 𝘅[2], "\t", 𝘅[3], "\t",
                  norm(𝘅[1:2]), "\t", λ, "\t", H, "\t", K);
    end
  end
  close(iod);

  ioe = open(path * "/error-t$(t_id).txt", "w");
  println(ioe, "Er\tEλ\tEH");
  Eλ = √Eλ; EH = √EH;
  println(ioe, Er, "\t", Eλ, "\t", EH);
  close(ioe);
  
  return Nothing;
end # analyze_output



