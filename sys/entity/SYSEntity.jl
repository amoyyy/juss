###############################################################################
## SYSEntity
###############################################################################
module SYSEntity
    using MaterialLib

    ###############################################################################
    abstract type Component end
    """
    Basic variables includes:
        name::Symbol
        keys::Dict{Symbol, DataType}
        data::Dict{DataType, Dict{Symbol, Any}}
                            # Dict{Symbol, Int64}
                            # Dict{Symbol, Symbol}
                            # Dict{Symbol, Float64}
                            # Dict{Symbol, Vector{Float64}}
    """
    #---------------------------------------------------------------------------
    # Several Universal Functions
    function Base.getindex(c::Component, tag::Symbol)
        if tag in keys(c.keys)
            return c.data[c.keys[tag]][tag]
        end
        nothing
    end

    function Base.getindex(c::Component, tag::String)
        return c[Symbol(tag)]
    end

    function add_value(c::Component, key::Symbol, data::Any)
        type = typeof(data)
        if key in keys(c.keys)
            println("ERROR: Variable: $tag => Value:$(c[key]) already exists!")
        else
            if !(type in keys(c.data))
                push!(c.data, typeof(data) => Dict{Symbol, type}())
            end
            push!(c.keys, (key => typeof(data)))
            push!(c.data[type], (key => data))
        end
    end
    #---------------------------------------------------------------------------
    include("components/element.jl")
    include("components/structure.jl")
    include("components/circuit.jl")
    ###############################################################################
    #
    abstract type Node end
    abstract type Face end

    #------------------------------------------------------------------------------
    function Base.getindex(p::Node, tag::Symbol)
        return p.data[tag]
    end

    function Base.getindex(p::Node, tag::String)
        return p.data[Symbol(tag)]
    end

    function Base.getindex(p::Face, tag::Symbol)
        return p.data[tag]
    end

    function Base.setindex!(p::Node, value::Float64, tag::Symbol)
        p.data[tag] = value
    end

    function Base.setindex!(p::Face, value::Float64, tag::Symbol)
        p.data[tag] = value
    end

    include("node/hydro_node.jl")
    include("node/thermal_node.jl")
    ###############################################################################
    #
    abstract type Mesh end
    include("meshing/hydro_1dmesh.jl")
    include("meshing/thermal_2dmesh.jl")
    ###############################################################################
    #
    include("system/plant.jl")
    ###############################################################################
end
