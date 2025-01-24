

@testset "Knot span" begin
  kv0 = KnotVector([0., 0., 0., 0.3, 0.5, 0.5, 0.6, 1., 1., 1.], 0, CLAMPED);
  @test get_knot_span_index(kv0, 0.0) == 3;
  @test get_knot_span_index(kv0, 0.5) == 6;
  @test get_knot_span_index(kv0, 1.0) == 7;

  kv1 = KnotVector([0., 0., 0., 1., 2., 3., 4., 4., 4.], 2, CLAMPED);
  @test get_knot_span_index(kv1, 0.0) == 3;
  @test get_knot_span_index(kv1, 0.7) == 3;
  @test get_knot_span_index(kv1, 1.0) == 4;
  @test get_knot_span_index(kv1, 1.3) == 4;
  @test get_knot_span_index(kv1, 2.0) == 5;
  @test get_knot_span_index(kv1, 2.1) == 5;
  @test get_knot_span_index(kv1, 3.0) == 6;
  @test get_knot_span_index(kv1, 3.9) == 6;
  @test get_knot_span_index(kv1, 4.0) == 6;
  @test_throws AssertionError get_knot_span_index(kv1, -0.1);
  @test_throws AssertionError get_knot_span_index(kv1,  4.1);

  kv1_dup = KnotVector([0., 0., 0., 1., 2., 3., 4., 4., 4.], 2, CLAMPED);
  @test kv1 == kv1_dup;
  @test hash(kv1) == hash(kv1_dup);

  kv1_nζ = KnotVector([0., 0., 0., 1.0001, 2., 3., 4., 4., 4.], 2, CLAMPED);
  @test kv1 != kv1_nζ;
  @test hash(kv1) != hash(kv1_nζ);

  @test_throws AssertionError KnotVector([2., 2., 2., 3., 4., 4., 4.], 3, CLAMPED);
  @test_throws AssertionError KnotVector([2., 2., 2., 3., 4., 4.], 2, CLAMPED);
  @test_throws AssertionError KnotVector([2., 2., 2., 3., 2.5, 4., 4., 4.], 2, CLAMPED);
  @test_throws AssertionError KnotVector([2., 2., 2., 3., 4., 4., 4.], 2, CLOSED);

  kv2 = KnotVector(4, 3, CLOSED);
  @test_throws AssertionError get_knot_span_index(kv2, -0.75);
  @test_throws AssertionError get_knot_span_index(kv2,  1.1 );
  @test_throws AssertionError get_knot_span_index(kv2,  1.75);
  @test get_knot_span_index(kv2,  0.0 ) ==  4;
  @test get_knot_span_index(kv2,  0.2 ) ==  4;
  @test get_knot_span_index(kv2,  0.4 ) ==  5;
  @test get_knot_span_index(kv2,  0.5 ) ==  6;
  @test get_knot_span_index(kv2,  0.75) ==  7;
  @test get_knot_span_index(kv2,  1.0 ) ==  7;

end

@testset "B-spline functions" begin
  # Ex 2.3 in Ref. [PT]
  kv1 = KnotVector([0., 0., 0., 1., 2., 3., 4., 4., 5., 5., 5.], 2, CLAMPED);
  @test norm(get_bspline_vals(kv1, 2.5) - [1/8, 6/8, 1/8]) < eps(8.);

  # from ALE fluids code
  kv2 = KnotVector([0., 0., 0., 1/3, 2/3, 1., 1., 1.], 2, CLAMPED);
  @test norm(get_bspline_vals(kv2, 0.0)       - [1., 0., 0.])   < eps(8.);
  @test norm(get_bspline_vals(kv2, 1/3-eps()) - [0., 1/2, 1/2]) < eps(8.);
  @test norm(get_bspline_vals(kv2, 1/3)       - [1/2, 1/2, 0.]) < eps(8.);

  @test kv2 == KnotVector(3, 2, CLAMPED);
  @test hash(kv2) == hash(KnotVector(3, 2, CLAMPED));

  # from https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-ex-1.html
  ζs3 = [0., 0., 0., 0.3, 0.5, 0.5, 0.6, 1., 1., 1.];
  kv3 = KnotVector(ζs3, 1, CLAMPED);
  @test norm(get_bspline_vals(kv3, 0.00) - [1.00, 0.00]) < eps(8.);
  @test norm(get_bspline_vals(kv3, 0.15) - [0.50, 0.50]) < eps(8.);
  @test norm(get_bspline_vals(kv3, 0.35) - [0.75, 0.25]) < eps(8.);
  @test norm(get_bspline_vals(kv3, 0.53) - [0.70, 0.30]) < eps(8.);
  @test norm(get_bspline_vals(kv3, 0.60) - [1.00, 0.00]) < eps(8.);
  @test norm(get_bspline_vals(kv3, 0.90) - [0.25, 0.75]) < eps(8.);
end

@testset "B-spline derivatives" begin
  # Ex 2.3 in Ref. [PT]
  kv1 = KnotVector([0., 0., 0., 1., 2., 3., 4., 4., 5., 5., 5.], 2, CLAMPED);
  @test norm(get_bspline_ders(kv1, 2.5, 2)[:,1] - [1/8, 6/8, 1/8]) < eps(8.);

  # from ALE fluids code
  kv2 = KnotVector([0., 0., 0., 1/3, 2/3, 1., 1., 1.], 2, CLAMPED);
  @test norm(get_bspline_ders(kv2, 0.0, 0)[:,1]   - [1., 0., 0.])    < eps(8.);
  @test norm(get_bspline_ders(kv2, 1/3-eps(), 0) .- [0., 1/2, 1/2])  < eps(8.);
  @test norm(get_bspline_ders(kv2, 1/3, 1)[:,1]   - [1/2, 1/2, 0.])  < eps(8.);
  @test norm(
    get_bspline_ders(kv2, 1/2, 1) - [0.125 -1.5; 0.75 0; 0.125 1.5]) < eps(8.);

  # from https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-ex-1.html
  ζs3 = [0., 0., 0., 0.3, 0.5, 0.5, 0.6, 1., 1., 1.];
  kv3 = KnotVector(ζs3, 1, CLAMPED);
  @test norm(get_bspline_ders(kv3, 0.00, 1)[:,2] - [-10/3, 10/3]) < eps(8.);
  @test norm(get_bspline_ders(kv3, 0.15, 1)[:,2] - [-10/3, 10/3]) < eps(8.);
  @test norm(get_bspline_ders(kv3, 0.35, 1)[:,2] - [-5   ,  5])   < eps(8.);
  @test norm(get_bspline_ders(kv3, 0.53, 1)[:,2] - [-10  , 10])   < eps(16.);
  @test norm(get_bspline_ders(kv3, 0.60, 1)[:,2] - [-2.5 , 2.5])  < eps(8.);
  @test norm(get_bspline_ders(kv3, 0.90, 1)[:,2] - [-2.5 , 2.5])  < eps(8.);

  kv4 = KnotVector(ζs3, 2, CLAMPED);
  @test norm(
    get_bspline_ders(kv4, 0.00, 2)[:,3] - [200/9, -320/9, 40/3])    < eps(64.);
  @test norm(
    get_bspline_ders(kv4, 0.15, 2)[:,3] - [200/9, -320/9, 40/3])    < eps(64.);
  @test norm(get_bspline_ders(kv4, 0.35, 2)[:,3] - [20, -70, 50])   < eps(64.);
  @test norm(get_bspline_ders(kv4, 0.53, 2)[:,3] - [200, -240, 40]) < eps(2e3);
  @test norm(
    get_bspline_ders(kv4, 0.60, 2)[:,3] - [10, -22.5, 12.5])        < eps(128.);
  @test norm(
    get_bspline_ders(kv4, 0.90, 2)[:,3] - [10, -22.5, 12.5])        < eps(128.);
  @test norm(
    get_bspline_ders(kv4, 1.00, 2)[:,3] - [10, -22.5, 12.5])        < eps(128.);

  @test_throws AssertionError get_bspline_ders(kv4, 0.50, -1);
  @test_throws AssertionError get_bspline_ders(kv4, 0.50,  3);
  @test_throws AssertionError get_bspline_ders(kv4,-0.01,  2);
  @test_throws AssertionError get_bspline_ders(kv4, 1.01,  2);
end

@testset "ζ to global indices" begin

  # from https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-ex-1.html
  ζs1 = [0., 0., 0., 0.3, 0.5, 0.5, 0.6, 1., 1., 1.];
  kv1 = KnotVector(ζs1, 2, CLAMPED);
  @test get_bspline_indices(kv1, 0.00) == [1, 2, 3];
  @test get_bspline_indices(kv1, 0.29) == [1, 2, 3];
  @test get_bspline_indices(kv1, 0.30) == [2, 3, 4];
  @test get_bspline_indices(kv1, 0.49) == [2, 3, 4];
  @test get_bspline_indices(kv1, 0.50) == [4, 5, 6];
  @test get_bspline_indices(kv1, 0.59) == [4, 5, 6];
  @test get_bspline_indices(kv1, 0.60) == [5, 6, 7];
  @test get_bspline_indices(kv1, 0.69) == [5, 6, 7];
  @test get_bspline_indices(kv1, 1.00) == [5, 6, 7];

  # from https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-ex-1.html
  # note basis functions number 1 and 8 are never returned
  kv2 = KnotVector(ζs1, 1, CLAMPED);
  @test get_bspline_indices(kv2, 0.00) == [2, 3];
  @test get_bspline_indices(kv2, 0.29) == [2, 3];
  @test get_bspline_indices(kv2, 0.30) == [3, 4];
  @test get_bspline_indices(kv2, 0.49) == [3, 4];
  @test get_bspline_indices(kv2, 0.50) == [5, 6];
  @test get_bspline_indices(kv2, 0.59) == [5, 6];
  @test get_bspline_indices(kv2, 0.60) == [6, 7];
  @test get_bspline_indices(kv2, 0.69) == [6, 7];
  @test get_bspline_indices(kv2, 1.00) == [6, 7];

end

@testset "Control point calculation" begin
  kvs = [KnotVector([0., 0., 0., 1/3, 2/3, 1., 1., 1.], 2, CLAMPED),
         KnotVector([0., 0., 0., 0., 1/3, 2/3, 1., 1., 1., 1.], 3, CLAMPED)];
  for kv in kvs
    for f in [x -> x, x -> 3*x+1, x -> 1-4*(x-0.5)^2]
      cps = get_1d_bspline_cps(kv, f);
      for ζ=0.0:0.1:1.0
        @test dot(get_bspline_vals(kv, ζ),
                  cps[get_bspline_indices(kv, ζ)]) ≈ f(ζ);
      end
    end
  end

  @test_throws AssertionError get_1d_bspline_cps(KnotVector([1., 1., 2., 3., 3.], 1, CLAMPED), x -> x);
  @test_throws AssertionError get_1d_bspline_cps(KnotVector([1., 1., 1., 1., 1., 2., 3., 3., 3., 3., 3.], 4, CLAMPED), x -> x);
  @test_throws AssertionError get_1d_bspline_cps(KnotVector([1., 1., 1., 2., 2., 3., 3., 3.], 2, CLAMPED), x -> x);
end

@testset "Unique 1d elements" begin
  kv1 = KnotVector(2, 2, CLAMPED);
  uel_num, num_el, uel_ids, uel_list = get_unique_1d_elements(kv1);
  @test uel_num  == 2; @test num_el == 2;
  @test uel_ids  == [1, 2];
  @test uel_list == [(0.0, 1/2), (1/2, 1.0)];

  kv2 = KnotVector(5, 2, CLAMPED);
  uel_num, num_el, uel_ids, uel_list = get_unique_1d_elements(kv2);
  @test uel_num  == 5; @test num_el == 5;
  @test uel_ids  == [1, 2, 3, 4, 5];
  @test uel_list == [(0.0, 1/5), (1/5, 2/5), (2/5, 3/5), (3/5, 4/5), (4/5, 1.0)];

  kv3 = KnotVector(8, 2, CLAMPED);
  uel_num, num_el, uel_ids, uel_list = get_unique_1d_elements(kv3);
  @test uel_num  == 5; @test num_el == 8;
  @test uel_ids  == [1, 2, 3, 3, 3, 3, 4, 5];
  @test uel_list == [(0.0, 1/8), (1/8, 2/8), (2/8, 3/8), (6/8, 7/8), (7/8, 1.0)];

  kv4 = KnotVector([0., 0., 0., 1., 2., 3., 4., 5., 6., 7., 8., 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9., 10., 11., 12., 13., 14., 15., 16., 17., 17., 17.], 2, CLAMPED);
  uel_num, num_el, uel_ids, uel_list = get_unique_1d_elements(kv4);
  @test uel_num  == 15; @test num_el == 26;
  @test uel_ids  == [1, 2, 3, 3, 3, 3, 4, 5, 6, 7, 8, 8, 8, 8, 8, 8, 9, 10, 11, 12, 13, 13, 13, 13, 14, 15];
  @test uel_list == [(0., 1.), (1., 2.), (2., 3.), (6., 7.), (7., 8.), (8., 8.1), (8.1, 8.2), (8.2, 8.3), (8.8, 8.9), (8.9, 9.), (9., 10.), (10., 11.), (11., 12.), (15., 16.), (16., 17.)];

  kv6 = KnotVector(3, 3, CLAMPED);
  uel_num, num_el, uel_ids, uel_list = get_unique_1d_elements(kv6);
  @test uel_num  == 3; @test num_el == 3;
  @test uel_ids  == [1, 2, 3];
  @test uel_list == [(0.0, 1/3), (1/3, 2/3), (2/3, 1.0)];

  kv7 = KnotVector(6, 3, CLAMPED);
  uel_num, num_el, uel_ids, uel_list = get_unique_1d_elements(kv7);
  @test uel_num  == 6; @test num_el == 6;
  @test uel_ids  == [1, 2, 3, 4, 5, 6];
  @test uel_list == [(0.0, 1/6), (1/6, 2/6), (2/6, 3/6), (3/6, 4/6), (4/6, 5/6), (5/6, 1.0)];

  kv8 = KnotVector(9, 3, CLAMPED);
  uel_num, num_el, uel_ids, uel_list = get_unique_1d_elements(kv8);
  @test uel_num  == 7; @test num_el == 9;
  @test uel_ids  == [1, 2, 3, 4, 4, 4, 5, 6, 7];
  @test uel_list == [(0.0, 1/9), (1/9, 2/9), (2/9, 3/9), (3/9, 4/9), (6/9, 7/9), (7/9, 8/9), (8/9, 1.0)];

  kv9 = KnotVector([0., 0., 0., 0., 1., 2., 3., 4., 5., 6., 7., 8., 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9., 10., 11., 12., 13., 14., 15., 16., 17., 17., 17., 17.], 3, CLAMPED);
  uel_num, num_el, uel_ids, uel_list = get_unique_1d_elements(kv9);
  @test uel_num  == 21; @test num_el == 26;
  @test uel_ids  == [1, 2, 3, 4, 4, 5, 6, 7, 8, 9, 10, 11, 11, 11, 11, 12, 13, 14, 15, 16, 17, 18, 18, 19, 20, 21];
  @test uel_list == [(0., 1.), (1., 2.), (2., 3.), (3., 4.), (5., 6.), (6., 7.), (7., 8.), (8., 8.1), (8.1, 8.2), (8.2, 8.3), (8.3, 8.4), (8.7, 8.8), (8.8, 8.9), (8.9, 9.), (9., 10.), (10., 11.), (11., 12.), (12., 13.), (14., 15.), (15., 16.), (16., 17.)];

  kv10 = KnotVector(5, 2, CLOSED);
  uel_num, num_el, uel_ids, uel_list = get_unique_1d_elements(kv10);
  @test uel_num  == 1;
  @test num_el   == 5;
  @test uel_ids  == [1, 1, 1, 1, 1];
  @test uel_list == [(0, 0.2)];

end

@testset "Collocate ζ" begin
  @test collocate_ζ(KnotVector(5, 2, CLOSED)) ≈ [0.1, 0.3, 0.5, 0.7, 0.9];
  @test collocate_ζ(KnotVector(5, 3, CLOSED)) ≈ [0.1, 0.3, 0.5, 0.7, 0.9];
end

