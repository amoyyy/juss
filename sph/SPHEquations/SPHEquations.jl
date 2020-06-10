module SPHEquations
    using Base.Threads
    using LinearAlgebra
    using SPHParticles
    using SPHKernels
    using MaterialLib

################################################################################
## Equation
################################################################################
    abstract type Equation end
    function initialize(T::Equation, P::SPHParticles.ParticlePool) end
    function loop(T::Equation, P::SPHParticles.ParticlePool, k::SPHKernels.Kernel) end
    function post_loop(T::Equation, P::SPHParticles.ParticlePool) end

    include("state.jl")
    include("continuity.jl")
    include("momentum.jl")
    include("energy.jl")
    include("heat_conduction.jl")
    include("boundary.jl")

################################################################################
## Group
################################################################################
    # DEBUG
    mutable struct Group
        ranks::Vector{Int64}
        dsts::Dict{Int64, Symbol}
        equations::Dict{Int64, Vector{Equation}}
        Group() = new(Vector{Int64}(), Dict{Int64, Symbol}(), Dict{Int64, Vector{Equation}}())
    end

    function add(g::Group, rank::Int64, eq::Equation)
        if rank in keys(g.ranks)
            if eq.dst == g.dsts[rank]
                push!(g.equations[rank], eq)
            else
                println("ERROR: Equation destination not match!!! $(eq)")
                exit()
            end
        else
            push!(g.ranks, rank)
            sort!(g.ranks)
            push!(g.equations, (rank => Vector{Equation}()))
            push!(g.equations[rank], eq)
            push!(g.dsts, rank=>eq.dst)
        end
    end

    function Base.display(g::Group)
        println("Equations:")
        for rank in g.ranks
            foreach(x->println("Rank: $rank , Destination: $(g.dsts[rank])   $x"), g.equations[rank])
        end
    end

    # in future ---> independent parallel handler
    function calculate_force(g::SPHEquations.Group, pools::SPHParticles.ParticlePoolUnion, kernel::SPHKernels.Kernel)
        # DEBUG solve equations
        #******************************************************************#
        for rank in g.ranks
            dst = g.dsts[rank]
            # initialize
            for equation in g.equations[rank]
                foreach(pool->initialize(equation, pool), pools[dst])
            end
            # loop
            nnps_a = SPHParticles.ParticlePool(:ALL)
            for pool in pools[dst]
                #@threads
                for d_idx in pool.pidx
                    #nnps_f = SPHParticles.ParticlePool(:Fluid)
                    #nnps_s = SPHParticles.ParticlePool(:Solid)
                    #SPHParticles.get_adjacency_particles(d_idx, pool, nnps_f)
                    #SPHParticles.get_adjacency_particles(d_idx, pool, nnps_s)
                    SPHParticles.get_adjacency_particles(d_idx, pool, nnps_a)
                    # DEBUG Better Interface Needed
                    for equation in g.equations[rank]
                        loop(equation, kernel, nnps_a, pool, d_idx)
                    end
"""                    for equation in g.equations[rank]
                        if equation.src == :NONE
                            loop(equation, pool, d_idx)
                        elseif equation.src == :Fluid
                            loop(equation, kernel, nnps_f, pool, d_idx)
                        elseif equation.src == :Solid
                            loop(equation, kernel, nnps_s, pool, d_idx)
                        elseif equation.src == :ALL
                            loop(equation, kernel, nnps_a, pool, d_idx)
                        end
                    end"""
                end
            end
            # post_loop
            for equation in g.equations[rank]
                foreach(pool->post_loop(equation, pool), pools[dst])
            end
        end
        #******************************************************************#
    end
################################################################################
end
