###############################################################################
## SYSSolver
###############################################################################
module SYSSolver
    using Base.Threads
    using SparseArrays
    using LinearAlgebra
    using KrylovMethods

    using SYSEntity
    using MaterialLib

###############################################################################
    const GRAVITY = 9.79
    const RANKS = Dict([(:HYDRO => 1),
                        (:THERMAL => 2)])
    # Numerical Models
    include("models/flow_resistence.jl")
    include("models/heat_transfer.jl")

###############################################################################
    abstract type BasicSolver end
    #
    include("hydro_solver/hydro_newton.jl")
    include("hydro_solver/hydro_jfnk.jl")
    include("thermal_solver/thermal_basic.jl")
    #

###############################################################################
# MainSolver
###############################################################################
    mutable struct Solver
        # solvers
        solver_ranks::Vector{Int64}
        sub_solvers::Dict{Int64, Vector{BasicSolver}}

        # Constructors
        Solver() = new(Vector{Int64}(), Dict{Int64, Vector{BasicSolver}}())
        Solver(p::SYSEntity.Plant, settings::NamedTuple) = init(Solver(), p, settings)
    end

    function init(s::Solver, p::SYSEntity.Plant, settings::NamedTuple)
        #------------------------------------------------------------------
        # HydroSolver
        for circuit in keys(p.circuits)
            add_solver(s, RANKS[:HYDRO], HydroJFNK(p, circuit))  # JFNK
            #add_solver(s, RANKS[:HYDRO], HydroNewton(p, circuit))  # Newton
        end
        #------------------------------------------------------------------
        # ThermalSolver
        for structure in keys(p.structures)
            add_solver(s, RANKS[:THERMAL], ThermalBasic(p, structure))
        end
        #------------------------------------------------------------------
        # Other Settings

        #------------------------------------------------------------------
        return s
    end

    function add_solver(s::Solver, rank::Int64, sub_s::BasicSolver)
        if rank in s.solver_ranks
            push!(s.sub_solvers[rank], sub_s)
        else
            push!(s.solver_ranks, rank)
            sort!(s.solver_ranks)
            push!(s.sub_solvers, (rank => Vector{BasicSolver}()))
            push!(s.sub_solvers[rank], sub_s)
        end
    end

    function set_dt(s::Solver, dt::Float64)
        for rank in s.solver_ranks
            foreach(solver->set_dt(solver, dt), s.sub_solvers[rank])
        end
    end

    # Optimilization to be fulfilled
    function solve(s::Solver)
        for rank in s.solver_ranks
            @threads for solver in s.sub_solvers[rank]
                @time solve(solver)
            end
        end
    end

    function solve(s::Solver, rank::Int64)
        for solver in s.sub_solvers[rank]
            @time solve(solver)
        end
    end

    function get_grids(s::Solver, structure::Symbol)
        for solver in s.sub_solvers[RANKS[:THERMAL]]
            if structure==solver.name
                return solver.mesh
            end
        end
        nothing
    end

    function get_cvs(s::Solver, ele::Symbol)
        for solver in s.sub_solvers[RANKS[:HYDRO]]
            if ele in keys(solver.mesh.circuit)
                return solver.mesh[ele]
            end
        end
        nothing
    end

"""###############################################################################
## Coupling Solver / temporarily removed...
###############################################################################
    # Node-Node Coupling
    struct NNCoupler{T1<:SYSEntity.Node, T2<:SYSEntity.Node} <: BasicSolver
        name::Symbol
        couplers::Vector{Pair{T1, T2}}
        #
        NNCoupler(T1::DataType, T2::DataType, name::Symbol=:NONE) =
        new{T1, T2}(name, Vector{Pair{T1, T2}}())
    end

    # Node-Nodes Coupling
    struct NNSCoupler{T1<:SYSEntity.Node, T2<:SYSEntity.Node} <: BasicSolver
        name::Symbol
        couplers::Vector{Pair{T1, Vector{T2}}}
        #
        NNSCoupler(T1::DataType, T2::DataType) =
        new{T1, T2}(:NONE, Vector{Pair{T1, Vector{T2}}}())
    end

    function add_coupler(s::NNCoupler, pairs)
        push!(s.couplers, pairs)
    end

    function add_coupler(s::NNSCoupler, pairs)
        push!(s.couplers, pairs)
    end"""
###############################################################################
end
