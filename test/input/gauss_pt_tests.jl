

@testset "Gauss points ξ" begin
  gpξ3 = GaussPointsξ(3);
  @test gpξ3.ngp == 3;
  @test gpξ3.ξs == [-sqrt(3/5), 0, sqrt(3/5)];
  @test gpξ3.ws == [5/9, 8/9, 5/9];

  gpξ4 = GaussPointsξ(4);
  @test gpξ4.ngp == 4;
  @test norm(gpξ4.ξs - [-sqrt(3/7+2/7*sqrt(6/5)), -sqrt(3/7-2/7*sqrt(6/5)), sqrt(3/7-2/7*sqrt(6/5)), sqrt(3/7+2/7*sqrt(6/5))]) < eps(1.);
  @test gpξ4.ws == [(18-sqrt(30))/36, (18+sqrt(30))/36, (18+sqrt(30))/36, (18-sqrt(30))/36];

  gpζ3 = GaussPointsζ(gpξ3, 0.5, 3.5);
  @test gpζ3.ngp == 3;
  @test norm(gpζ3.ζs - (gpξ3.ξs * 1.5 .+  2.0)) < eps(1.);
  @test norm(gpζ3.ws - gpξ3.ws * 1.5) < eps(1.);

  gpζ4 = GaussPointsζ(gpξ4, 1.0, 6.0);
  @test gpζ4.ngp == 4;
  @test norm(gpζ4.ζs - (gpξ4.ξs * 2.5 .+  3.5)) < eps(1.);
  @test norm(gpζ4.ws - gpξ4.ws * 2.5) < eps(1.);
end

