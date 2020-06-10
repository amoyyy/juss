###############################################################################
## Solver
###############################################################################
module SPHSolver
    using SPHIntegrators
    using SPHKernels
    using SPHParticles
    using SPHEquations

    include("fluid_condition.judge.jl")

    mutable struct Solver{T1<:SPHIntegrators.Integrator, T2<:SPHKernels.Kernel}
        dim::Int
        #
        integrator::T1
        kernel::T2
        #
        tf::Float64
        tc::Float64
        dt::Float64
        # the ratio of smoothing length and delta_r
        hdx::Float64
        #
        opt_freq::Int
        output_at_times::Vector{Float64}
        output_dir::String
        #
        particles::SPHParticles.ParticlePoolUnion
        equations::SPHEquations.Group

        Solver(dim::Int64, integrator::SPHIntegrators.Integrator, kernel::SPHKernels.Kernel, tf::Float64=1.0, tc::Float64=0.0, dt::Float64=1e-3, hdx::Float64=3.1, freq::Int64=100, output_at_times::Vector{Float64}=[1e-1, 1.0]) =
        new{typeof(integrator), typeof(kernel)}(dim, integrator, kernel, tf, tc, dt, hdx, freq, output_at_times)

        #Solver(dim, integrator, kernel) = Solver(dim, integrator, kernel, 1.0, 0.0, 1e-3, 3.1, 100, [1e-1, 1.0])
    end

################################################################################
    # DEBUG
    function pre_solve(sol::Solver, pools::SPHParticles.ParticlePoolUnion)
    end

    # DEBUG
    function solve(sol::Solver)
        iteration_times::Int = 0
        println("t=$(sol.tc)")
        output(sol)
        while(sol.tc < sol.tf)
            iteration_times += 1
            @time SPHIntegrators.timestep(sol.integrator, sol.dt)
            @time update(sol)
            @time output(sol, iteration_times)
            println("t=$(sol.tc)")
        end
    end


    # exec particles transfer and update the solver settings
    function update(s::Solver)
        print("     update:         ")
        #******************************************************************#
        # transfer particles : pools
        SPHParticles.update(s.particles)
        # update time
        s.tc += s.dt
        # update time_step
        #if s.adaptive_timestep
        #    #s.dt = FluidConditionJudge.cal_timestep(s.dt, s.hdx, pools["Fluid"])
        #end
        # update smooth length
        #if !s.fixed_h
        #    s.hdx = cal_smoothing_length(s.dim, s.hdx, pools["Fluid"])
        #end
        #******************************************************************#
    end

################################################################################
    function set_target(sol::Solver, eqs::SPHEquations.Group, pools::SPHParticles.ParticlePoolUnion)
        sol.particles = pools
        sol.equations = eqs
        sol.integrator.eqs = eqs
        sol.integrator.kernel = sol.kernel
        sol.integrator.particles = sol.particles
    end


    function set_output_dir(sol::Solver, path::String)
        if isdir(path)
            rm(path, recursive=true)
        end
        mkpath(path)
        sol.output_dir = path
    end

################################################################################
    function output(sol::Solver, iteration_times::Int=0)
        print("     output:         ")
        if mod(iteration_times, sol.opt_freq)==0 || sol.tc in sol.output_at_times
            dirpath = string("$(sol.output_dir)/t=$(sol.tc)s")
            for pools in sol.particles[]
                foreach(x->SPHParticles.output(x, dirpath, "$(x.tag)_$(x.id).nc"), pools)
            end
        end
    end

################################################################################
end
