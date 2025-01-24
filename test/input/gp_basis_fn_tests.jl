


@testset "GP basis functions ζ" begin
  kv1  = KnotVector(4, 2, CLAMPED);
  fn1ζ = GpBasisFnsζ(0.1, 0.4, kv1);
  @test fn1ζ.w   == 0.1;
  @test fn1ζ.N   == get_bspline_ders(kv1, 0.4, 2)[:,1];
  @test fn1ζ.dN  == get_bspline_ders(kv1, 0.4, 2)[:,2];
  @test fn1ζ.ddN == get_bspline_ders(kv1, 0.4, 2)[:,3];

  @test_throws AssertionError GpBasisFnsζ(0.1, -0.1, kv1);
  @test_throws AssertionError GpBasisFnsζ(0.1,  1.1, kv1);
end


@testset "GP basis functions ζα" begin

  kv3  = KnotVector([0., 0., 0., 1., 3., 3.5, 6., 6., 6.], 2, CLAMPED);
  kv4  = KnotVector([2., 2., 2., 4., 5., 5., 5.], 2, CLAMPED);

  for ζ1=0.0:1.5:6.0, ζ2=2.0:1.5:5.0
    w1 = rand(); w2 = rand();
    fn1 = GpBasisFnsζ(w1, ζ1, kv3);
    fn2 = GpBasisFnsζ(w2, ζ2, kv4);
    fnα = GpBasisFnsζα(fn1, fn2);
    @test fnα.w == fn1.w * fn2.w;
    for id1=1:3, id2=1:3
      @test fnα.N[id1 + 3*(id2-1)]        == fn1.N[id1]   * fn2.N[id2];
      @test fnα.∂Nα[id1 + 3*(id2-1), 1]   == fn1.dN[id1]  * fn2.N[id2];
      @test fnα.∂Nα[id1 + 3*(id2-1), 2]   == fn1.N[id1]   * fn2.dN[id2];
      @test fnα.∂∂Nαβ[id1 + 3*(id2-1), 1] == fn1.ddN[id1] * fn2.N[id2];
      @test fnα.∂∂Nαβ[id1 + 3*(id2-1), 2] == fn1.N[id1]   * fn2.ddN[id2];
      @test fnα.∂∂Nαβ[id1 + 3*(id2-1), 3] == fn1.dN[id1]  * fn2.dN[id2];
    end
  end
end


@testset "Boundary GP basis functions" begin

  kv3_short  = KnotVector([1., 1., 1., 2., 3., 4., 5., 6., 7., 8., 8., 8.], 2, CLAMPED);
  kv3_long   = KnotVector([1., 1., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 12., 12.], 2, CLAMPED);
  fns3_short = LineGpBasisFns(kv3_short, 3);
  fns3_long  = LineGpBasisFns(kv3_long,  3);

  @test fns3_short.nel ==  7;
  @test fns3_long.nel  == 11;
  @test maximum(fns3_short.uel_ids) == 5;
  @test maximum(fns3_long.uel_ids)  == 5;

  for i_fn=1:5, j_gp=1:3
    @test fns3_short.ufns[i_fn, j_gp] == fns3_long.ufns[i_fn, j_gp];
  end

  ## manually check some of the basis functions and derivatives
  gpζ01 = GaussPointsζ(GaussPointsξ(3), 0.0, 1.0);
  @test fns3_short.ufns[1,1].N   == get_bspline_ders(kv3_short, gpζ01.ζs[1]+1, 0)[:,1];
  @test fns3_short.ufns[2,2].dN  == get_bspline_ders(kv3_long,  gpζ01.ζs[2]+2, 1)[:,2];
  @test fns3_short.ufns[4,3].N   == get_bspline_ders(kv3_long,  gpζ01.ζs[3]+8, 0)[:,1];
  @test fns3_short.ufns[4,3].ddN == get_bspline_ders(kv3_long,  gpζ01.ζs[3]+8, 2)[:,3];
  @test fns3_short.ufns[5,1].N   == get_bspline_ders(kv3_long,  gpζ01.ζs[1]+11,0)[:,1];

  # asymmetric knot vector
  kv2  = KnotVector([0., 0., 0., 3., 4., 4.5, 6.5, 9., 9.6, 10.1, 10.1, 10.1], 2, CLAMPED);
  fns2 = LineGpBasisFns(kv2, 3);
  @test maximum(fns2.uel_ids) == 7;
  @test fns2.ufns[1,1].N != fns2.ufns[7,3].N;
  @test fns2.ufns[1,1].N != fns2.ufns[7,1].N;
  gpζ03 = GaussPointsζ(GaussPointsξ(3), 0.0, 3.0);
  @test fns2.ufns[1,1].N == get_bspline_ders(kv2, gpζ03.ζs[1], 0)[:,1];
  gpζ96 = GaussPointsζ(GaussPointsξ(3), 9.0, 9.6);
  @test fns2.ufns[6,2].N == get_bspline_ders(kv2, gpζ96.ζs[2], 0)[:,1];
end


