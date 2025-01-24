

@testset "Geo dyn stress" begin

  @testset "Curved patch" begin
    # currently, to avoid namespace issues with unit tests, we define the
    # Params struct twice
    p2 = MembraneAleFem.Params(
           motion=MembraneAleFem.STATIC,
           scenario=MembraneAleFem.F_CAVI,
           num1el=2, num2el=3, output=false);

    mesh = Mesh(Params(motion=STATIC, scenario=F_CAVI,
                       num1el=2, num2el=3, output=false));

    xms  = zeros(Float64, mesh.numnp, XDIM);      # positions
    cps  = zeros(Float64, mesh.numnp, mesh.ndf);  # control points

    # initialize position
    A = 2.2; B = 1.7;
    x_cps = get_1d_bspline_cps(mesh.kv1, ζ1 -> 4*ζ1+2);
    y_cps = get_1d_bspline_cps(mesh.kv2, ζ2 -> 2*ζ2-3);
    z_cps = get_2d_bspline_cps(mesh.kv1, mesh.kv2,
                               (ζ1, ζ2) -> A*(4*ζ1-2)^2 + B*(2*ζ2-1)^2);

    for node=1:mesh.numnp
      xms[node, Int64(Dof.xm)] = x_cps[(node-1) % mesh.num1np + 1];
      xms[node, Int64(Dof.ym)] = y_cps[floor(Int64, (node-1)/mesh.num1np) + 1];
      xms[node, Int64(Dof.zm)] = z_cps[node];
    end

    # loop over elements
    for el_id=1:mesh.numel
      xms_el = complex(xms[mesh.IX[:,el_id], :]);
      cps_el = complex(cps[mesh.IX[:,el_id], :]);

      # loop over Gauss points
      for gp_id = 1:GP1D^ζDIM
        N     = get_N(el_id, gp_id, mesh);
        ∂Nα   = get_∂Nα(el_id, gp_id, mesh);
        ∂∂Nαβ = get_∂∂Nαβ(el_id, gp_id, mesh);
        gpw   = get_gpw(el_id, gp_id, mesh);
        gds   = GeoDynStress(xms_el, cps_el, mesh.dofs, N, ∂Nα, ∂∂Nαβ, p2);

        x = gds.𝘅[1]; y = gds.𝘅[2]; z = gds.𝘅[3];
        @test z ≈ A * (x - 4)^2 + B * (y + 2)^2;
        @test norm(gds.𝗮_α .- [4. 0.; 0. 2.; 8*A*(x-4) 4*B*(y+2)]) < eps(4.e2);
        @test norm(gds.a△αβ * (transpose(gds.𝗮_α) * gds.𝗮_α) - [1. 0.; 0. 1.]) < eps(1.e2)
        @test norm((transpose(gds.𝗮_α) * gds.𝗮_α) * gds.a△αβ - [1. 0.; 0. 1.]) < eps(1.e2)
        @test abs( tr(gds.a△αβ * gds.b_αβ)/2 - gds.H) < eps(1.e1);
        @test abs(det(gds.a△αβ * gds.b_αβ)   - gds.K) < eps(1.e1);
      end # Gauss point loop
    end # element loop

  end # curved patch

  @testset "Flat patch" begin
    # currently, to avoid namespace issues with unit tests, we define the
    # Params struct twice
    p2 = MembraneAleFem.Params(
           motion=MembraneAleFem.STATIC,
           scenario=MembraneAleFem.F_CAVI,
           num1el=3, num2el=2, output=false);

    mesh = Mesh(Params(motion=STATIC, scenario=F_CAVI,
                       num1el=3, num2el=2, output=false));

    xms   = zeros(Float64, mesh.numnp, XDIM);      # positions
    cps_c = zeros(Float64, mesh.numnp, mesh.ndf);  # Couette control points
    cps_p = zeros(Float64, mesh.numnp, mesh.ndf);  # Poiseuille control points

    # initialize position
    x_cps = get_1d_bspline_cps(mesh.kv1, ζ1 -> 4*ζ1+2);
    y_cps = get_1d_bspline_cps(mesh.kv2, ζ2 -> 2*ζ2-3);

    # initialize flat couette and poiseuille velocities
    ω = 12.;
    U = 16.;
    c_vs = get_1d_bspline_cps(mesh.kv2, ζ2 ->  2*ω*ζ2);
    p_vs = get_1d_bspline_cps(mesh.kv1, ζ1 -> -4*U*(ζ1^2 - ζ1));

    for node=1:mesh.numnp
      xms[node, Int64(Dof.xm)] = x_cps[(node-1) % mesh.num1np + 1];
      xms[node, Int64(Dof.ym)] = y_cps[floor(Int64, (node-1)/mesh.num1np) + 1];

      cps_c[node, Int64(Dof.vx)] = c_vs[floor(Int64, (node-1)/mesh.num1np) + 1];
      cps_p[node, Int64(Dof.vy)] = p_vs[(node-1) % mesh.num1np + 1];
    end

    area = 0.0;

    # loop over elements
    for el_id=1:mesh.numel
      xms_el  =   complex(xms[mesh.IX[:,el_id], :]);
      cps_cel = complex(cps_c[mesh.IX[:,el_id], :]);
      cps_pel = complex(cps_p[mesh.IX[:,el_id], :]);

      # loop over Gauss points
      for gp_id = 1:GP1D^ζDIM
        N     = get_N(el_id, gp_id, mesh);
        ∂Nα   = get_∂Nα(el_id, gp_id, mesh);
        ∂∂Nαβ = get_∂∂Nαβ(el_id, gp_id, mesh);
        gpw   = get_gpw(el_id, gp_id, mesh);
        gds_c = GeoDynStress(xms_el, cps_cel, mesh.dofs, N, ∂Nα, ∂∂Nαβ, p2);
        gds_p = GeoDynStress(xms_el, cps_pel, mesh.dofs, N, ∂Nα, ∂∂Nαβ, p2);

        area += gds_p.JΩ * gpw;

        # since gds_c and gds_p are generated from the same control point
        # positions, the position (x, y, z) is the same here
        x = gds_c.𝘅[1]; y = gds_c.𝘅[2]; z = gds_c.𝘅[3];

        @test norm(gds_c.𝞂 - p2.ζv * [0.; 0.; 2 * ω / 8]) < eps(2.e2);
        @test norm(gds_p.𝞂 - p2.ζv * [0.; 0.; U*(4-x)/8]) < eps(3.e2);
      end # Gauss point loop
    end # element loop
    @test area ≈ 8.0;

  end # flat patch

end

