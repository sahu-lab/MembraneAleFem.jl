


export Params

include("Enums.jl")


"""
    Params

Parameters that must be specified for any simulation.

For [`Scenario`](@ref)-specific data that must be included via keyword
arguments, see [`check_params`](@ref MembraneAleFem.check_params).
The following information is contained:
- `motion::`[`Motion`](@ref)     → type of motion
- `scenario::`[`Scenario`](@ref) → type of scenario
- `num1el::Int64`   → number of elements in ``ζ^1``
- `num2el::Int64`   → number of elements in ``ζ^2``
- `output::Bool`    → write output (yes by default)
- `length::Float64` → length, in real space
- `kb::Float64`     → mean bending modulus
- `kg::Float64`     → Gaussian bending modulus
- `ζv::Float64`     → membrane viscosity
- `pn::Float64`     → normal pressure (body force)
- `poly::Int64`     → polynomial order of basis functions
- `gp1d::Int64`     → number of Gauss points, in 1-D
- `nders::Int64`    → number of 1-D B-spline derivatives
- `nen::Int64`      → number of element nodes
- `αdb::Float64`    → Dohrmann--Bochev parameter
- `αm::Float64`     → mesh parameter
- `εk::Float64`     → tolerance numerical calculation of matrix K
- `εnr::Float64`    → Netwon--Raphson tolerance
"""
@with_kw struct Params
  motion::Motion       = ALEVB           # motion
  scenario::Scenario   = F_PULL          # scenario
  num1el::Int64        = 17              # number of elements in ζ1
  num2el::Int64        = 17              # number of elements in ζ2
  output::Bool         = true            # write output
  length::Float64      = 64.0            # length, in real space
  kb::Float64          =  1.0            # mean bending modulus
  kg::Float64          = -0.5            # Gaussian bending modulus
  ζv::Float64          =  1.0            # membrane viscosity
  pn::Float64          =  0.0            # normal pressure (body force)
  poly::Int64          = 2               # polynomial order of basis functions
  gp1d::Int64          = 3               # number of Gauss points, in 1-D
  nders::Int64         = 2               # number of 1-D B-spline derivatives
  nen::Int64           = (poly+1)^ζDIM   # number of element nodes
  αdb::Float64         = length^2        # Dohrmann--Bochev parameter
  αm::Float64          = 1.0             # mesh parameter
  εk::Float64          = 1.e-15;         # numerical calculation of matrix K
  εnr::Float64         = 1.e-12;         # Netwon--Raphson tolerance
end # Params



"""
    check_params(p::Params; args...)

Check relevant parameters, including scenario-specific keyword arguments.

For any simulation, the fields `t0`, `t0_id`, and `Δts` are required.
If data is output (as is the default, for `p.output`), the output path
`out_path` and file name `out_file` are necessary.
Finally, the information needed to solve the various [`Scenario`](@ref)s is
checked.
"""
function check_params(
    p::Params;
    args...
  )

  BLAS.set_num_threads(1);

  ks = keys(args);

  @assert :Δts   ∈ ks "\nneed list of time steps 'Δts'";
  @assert :t0    ∈ ks "\nneed initial time 't0'";
  @assert :t0_id ∈ ks "\nneed initial time ID 't0_id'";

  if p.output
    @assert :out_file ∈ ks "\nneed output file 'out_file'";
    @assert :out_path ∈ ks "\nneed output path 'out_path'";
    mkpath(args[:out_path]);
  end

  @assert p.nen == (p.poly+1)^ζDIM "\nincorrect number of element nodes";
  
  # the following parameters are required to use SVector and SMatrix, and so are
  # kept as global constants
  @assert p.nen   == NEN;
  @assert p.poly  == POLY;
  @assert p.nders == NDERS;

  if p.scenario == F_BEND
    @assert :bend_tm ∈ ks "\nneed final ramp-up time 'bend_tm'";
    @assert :bend_mf ∈ ks "\nneed final applied bending moment 'bend_mf'";
    @assert args[:bend_mf] == p.kb / 2 / p.length "\nα = 1/2 is the half-angle";
    @assert p.motion != STATIC "\nmesh cannot be static for $(F_BEND) scenario";
    @assert p.motion != ALEVB "\n$(F_BEND) not implemented for $(ALEVB) motion";
    @assert p.pn == 0.0 "$(F_BEND) with normal pressure not implemented";

  elseif p.scenario == F_PULL
    @assert :pull_speed ∈ ks "\nneed pull speed 'pull_speed'";

  elseif p.scenario == F_CAVI
    @assert p.motion == STATIC "\nneed static motion for cavity scenario";

  elseif p.scenario == F_COUE
    @assert p.motion == STATIC "\nneed static motion for Couette scenario";

  elseif p.scenario == F_POIS
    @assert p.motion == STATIC "\nneed static motion for Poiseuille scenario";

  end

  implemented_scenarios = [F_CAVI, F_COUE, F_POIS, F_PULL, F_BEND];
  if p.scenario ∉ implemented_scenarios
    @assert false "\n$(p.scenario) scenario not implemented!";
  end

  return Nothing;
end # check_params



"""
    display_params(p::Params; args...)

Display and serialize user-specified parameters and arguments.
"""
function display_params(
    p::Params;
    args...
  )

  if !p.output
    return Nothing;
  end

  open(args[:out_path] * "/params.txt", "a") do io
    show(io, p);
    println(io, "\n");
    show(io, args);
    println(io, "\n");
    println(io, "num BLAS  threads = ", BLAS.get_num_threads());
    println(io, "num Julia threads = ", Threads.nthreads());
  end

  serialize(args[:out_path] * "/params.dat",  p);
  serialize(args[:out_path] * "/args.dat", args);

  return Nothing;
end # display_params



# parameters
const POLY  = 2;  # basis function polynomial order
const GP1D  = 3;  # number of Gaussian quadrature points in 1-D
const NDERS = 2;  # number of basis function derivatives to calculate


## *** the following parameters should NEVER be changed *** ##

# dimensionality
const ζDIM = 2;
const XDIM = 3;

# number of element nodes
const NEN = (POLY+1)^ζDIM;

# structure of compressed tensors
const VOIGT = 3;

# number of element Dohrmann--Bochev degrees of freedom
const NEDBDF = 3;
