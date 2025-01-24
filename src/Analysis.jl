
include("analysis/GeoDynStress.jl");
include("analysis/FiniteElement.jl");
include("analysis/PullForce.jl");


"""
    run_analysis!(mesh, xms, cps, p::Params; args...)

Execute the finite element analysis and step forward in time through `args[:Δts]`.

At every time step, carry out Newton--Raphson iteration to advance to the next
time step.
The global residual vector and tangent diffusion matrix are generated at every
iteration, and mesh positions are updated according to the mesh velocities.
"""
function run_analysis!(
    mesh::Mesh,
    xms::Matrix{Float64},
    cps::Matrix{Float64},
    p::Params;
    args...
  )

  if p.output
    open(args[:out_path] * "/" * args[:out_file], "a") do io
      println(io, "\nStarting to run analysis...");
    end
  end

  times = [sum(args[:Δts][1:i]) for i=1:length(args[:Δts])] .+ args[:t0];

  # see if new output files are generated
  append = (:in_path ∈ keys(args)) && (args[:out_path] == args[:in_path]);

  # record simulation times, with index
  if p.output
    open(args[:out_path] * "/times.txt", "a") do io_time
      println(io_time, "$(length(args[:Δts])+args[:t0_id]) times:");
      for (t_id, time) in enumerate(times)
        println(io_time, "$(t_id+args[:t0_id])\t$(time)");
      end
    end
  end

  # pull force local mappings
  if p.scenario == F_PULL
    @assert length(get_v_order(mesh.dofs)) == XDIM "f_pull needs 3-D velocity";
    adj_el_ids, adj_node_map = get_adj_maps(mesh.num1el, mesh.numel,
                                            mesh.IX, p.poly);
    if !append && p.output
      open(args[:out_path] * "/f-pull.txt", "a") do io_pull
        println(io_pull, "time\txp\typ\tzp\tfx\tfy\tfz");
      end
    end
  end


  # time step loop
  for (t_id, Δt) in enumerate(args[:Δts])

    if p.output
      open(args[:out_path] * "/" * args[:out_file], "a") do io
        println(io, "\n-> Time ", times[t_id]);
      end
    end

    # update mesh positions: mesh velocity at time t is the initial guess for
    # the mesh velocity at time t + Δt
    update_xms!(p.motion, xms, cps, Δt, mesh.dofs);

    # step forward in time
    time_step!(mesh, xms, cps, times[t_id], Δt, p; args...);

    if p.output
      writedlm(args[:out_path] * "/t$(t_id+args[:t0_id])-xms.txt", xms);
      writedlm(args[:out_path] * "/t$(t_id+args[:t0_id])-cps.txt", cps);
    end

    # pull force calculation
    if p.scenario == F_PULL && p.output
      # assume all nodes from the pull element are moved the same amount
      x_pull   = sum(xms[mesh.IX[:, get_pull_el_id(mesh.numel)],:], dims=1);
      x_pull ./= size(mesh.IX)[1];
      f_pull   = calc_pull_force(mesh, xms, cps, adj_el_ids, adj_node_map, p);
      open(args[:out_path] * "/f-pull.txt", "a") do io_pull
        println(io_pull, times[t_id], "\t",
                         x_pull[1], "\t", x_pull[2], "\t", x_pull[3], "\t",
                         f_pull[1], "\t", f_pull[2], "\t", f_pull[3]);
      end
    end

  end # time step loop

  if p.output
    open(args[:out_path] * "/" * args[:out_file], "a") do io
      println(io, "\nCompleted running analysis...\n");
    end
  end

  return Nothing;
end # run_analysis


