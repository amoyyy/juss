    ###############################################################################
    # thermal
    ###############################################################################
"""    const BOUNDARY_TYPE = Dict([(:I => 1),
                                (:II => 2),
                                (:III => 3),
                                (:IV => 4)])"""

    mutable struct Grid <: Node
        id::Int64
        structure::Structure
        id_axial::Int64
        id_radial::Int64
        type::Int64
        #
        data::Dict{Symbol,Float64}
        links::Vector{Tuple{Grid,Float64,Float64}}
        boundary::CV
        #
        Grid(id::Int64, structure::Structure, axial::Int64, radial::Int64) =
        new(id, structure, axial, radial, 0, Dict{Symbol,Float64}(),
                 Vector{Tuple{Grid,Float64,Float64}}(), CV())
    end

    function init(g::Grid, mat::MaterialLib.Material)
        # geometry
        push!(g.data, (:cordin1 => (g.id_axial-0.5)*g.structure[:DELTA_X]))
        push!(g.data, (:delta_cordin1 => g.structure[:DELTA_X]))
        push!(g.data, (:cordin2 => (g.id_radial-0.5)*g.structure[:DELTA_R]+g.structure[:INNER_CORDIN]))
        push!(g.data, (:delta_cordin2 => g.structure[:DELTA_R]))
        push!(g.data, (:delta_v => g[:cordin2]*g[:delta_cordin1]*g[:delta_cordin2]))
        # property change with TEMPERATURE
        push!(g.data, (:t => g.structure[:TEMPERATURE][g.id]))
        push!(g.data, (:rho => MaterialLib.ro(mat, g[:t], 0.0)))
        push!(g.data, (:lamda => MaterialLib.lamda(mat, g[:t], 0.0)))
        push!(g.data, (:cp => MaterialLib.c_p(mat, g[:t], 0.0)))
        push!(g.data, (:heat_resistence1 => g[:delta_cordin1] / g[:lamda]))
        push!(g.data, (:heat_resistence2 => g[:delta_cordin2] / g[:lamda]))
        # other vars
        push!(g.data, (:q_sc => 0.0))
        push!(g.data, (:q_sp => 0.0))
        set_boundary(g)
    end

    function set_boundary(g::Grid)
        # boundary
        if g.id_radial == 1
            g.type = g.structure[:INNER_BOUNDARY_TYPE]
            push!(g.data, (:boundary_cordin => g.structure[:INNER_CORDIN]))
        elseif g.id_radial == size(g.structure)[2]
            g.type = g.structure[:OUTER_BOUNDARY_TYPE]
            push!(g.data, (:boundary_cordin => g.structure[:OUTER_CORDIN]))
        else
            g.type = 0
        end
        #
        if g.type == 1
            push!(g.data, (:twall => g.structure[:INNER_BOUNDARY_VALUE][g.id_axial]))
        elseif g.type == 2
            push!(g.data, (:qwall => g.structure[:INNER_BOUNDARY_VALUE][g.id_axial]))
        elseif g.type == 3 || g.type == 4
            push!(g.data, (:tfluid => 0.0))
            push!(g.data, (:h_coeff => 0.0))
        end
    end

    function fetch_boundary(g::Grid)
        # from cv to grid
        if g.type == 4
            g[:h_coeff] = g.boundary[:h_coeff]
            g[:tfluid] = g.boundary[:t]
            #g.boundary[:q] = g.boundary[:q_bak]
        end
    end

    function feedback_boundary(g::Grid)
        # from grid to cv
        if g.type == 4
            h_modify = 1.0/(1.0/g[:h_coeff] + 0.5*g[:heat_resistence2])
            g.boundary[:q] += h_modify*(g[:t]-g.boundary[:t])*2.0*pi*g[:boundary_cordin]/g.boundary.ele[:AREA]
            #println(g[:t],"  ",g.boundary[:t],"  ",h_modify)
        end
    end

    function link(g1::Grid, g2::Grid)
        dimension_coeff::Float64 = 0.0
        direction_coeff::Float64 = -1.0
        if g2.id_axial == g1.id_axial
            dimension_coeff = 1.0
        end
        if g2.id_axial > g1.id_axial || g2.id_radial > g1.id_radial
            direction_coeff = 1.0
        end
        push!(g1.links, (g2, dimension_coeff, direction_coeff))
        push!(g2.links, (g1, dimension_coeff, -direction_coeff))
    end

    function link(g1::Grid, g2::Nothing)
        nothing
    end

    function link(g1::Grid, g2::CV)
        g1.boundary = g2
        fetch_boundary(g1)
        feedback_boundary(g1)
    end
