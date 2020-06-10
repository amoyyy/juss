"""module SPHEquations
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
    include("heat_conduction.jl")
    include("boundary.jl")

################################################################################
## Group
################################################################################
    # DEBUG
    mutable struct Group
        ranks::Vector{Int64}
        equations::Dict{Int64, Vector{Equation}}
        Group() = new(Vector{Int64}(), Dict{Int64, Vector{Equation}}())
    end

    function add(g::Group, rank::Int64, eq::Equation)
        if rank in g.ranks
            push!(g.equations[rank], eq)
        else
            push!(g.ranks, rank)
            sort!(g.ranks)
            push!(g.equations, (rank => Vector{Equation}()))
            push!(g.equations[rank], eq)
        end
    end

    function Base.display(g::Group)
        println("Equations:")
        for rank in g.ranks
            foreach(x->println("Rank: $rank    $x"), g.equations[rank])
        end
    end

    function calculate_force(g::SPHEquations.Group, pools::SPHParticles.ParticlePoolUnion, kernel::SPHKernels.Kernel)
        # DEBUG solve equations
        for rank in g.ranks
            #******************************************************************#
            # initialize
            for equation in g.equations[rank]
                foreach(pool->initialize(equation, pool), pools[equation.dst])
            end
            # loop
            for equation in g.equations[rank]
                if equation.src == :NONE
                    @threads for d_idx in pools[equation.dst].pidx
                        loop(eq, pool, d_idx)
                    end
                else

                foreach(pool->loop_one_equation(equation, pool, kernel), pools[equation.dst])
            #end
            # post_loop
            #for equation in g.equations[rank]
                foreach(pool->post_loop(equation, pool), pools[equation.dst])
            end
            #******************************************************************#
        end
    end

    # in future ---> independent parallel handler
    function loop_one_equation(eq::SPHEquations.Equation, pool::SPHParticles.ParticlePool, kernel::SPHKernels.Kernel)
        if eq.src == :NONE
            @threads for d_idx in pool.pidx
                loop(eq, pool, d_idx)
            end
        else
            @threads for d_idx in pool.pidx
                nnps = SPHParticles.ParticlePool(eq.src)
                SPHParticles.get_adjacency_particles(d_idx, pool, nnps)
#                println(typeof(eq), "  ",d_idx," ",nnps.num_of_particles)
                if nnps.num_of_particles > 0
                    loop(eq, kernel, nnps, pool, d_idx)
                end
            end
        end
    end
################################################################################
end
"""
