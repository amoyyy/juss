    ##------------------------------------------------------------------------------
    ## Circuit
    ##------------------------------------------------------------------------------
    mutable struct Circuit
        name::Symbol
        fluid::Symbol
        keys::Vector{Symbol}
        elements::Dict{Symbol, Element}

        Circuit(name::Symbol) = new(name, :Sodium, Vector{Symbol}(), Dict{Symbol,Element}())
        Circuit(name::Symbol, input::NamedTuple) = init(Circuit(name), input)
    end

    function init(c::Circuit, input::NamedTuple)
        c.fluid = Symbol(input[:FLUID_MEDIA])
        for key in keys(input[:ELEMENT])
            if key in keys(c.elements)
                println("ERROR: Element: $(string(key)), already exists!")
            else
                add_element(c, Element(key, input[:ELEMENT][key]))
            end
        end
        for key in keys(input[:LINKTABLE])
            foreach(ele->add_next(c[key], c[ele]), input[:LINKTABLE][key])
        end
        return c
    end

    function add_element(c::Circuit, ele::Element)
        push!(c.keys, ele.name)
        push!(c.elements, (ele.name => ele))
    end

    #---------------------------------------------------------------------------------
    function Base.getindex(c::Circuit, ele::Symbol)
        return c.elements[ele]
    end

    function Base.getindex(c::Circuit, ele::String)
        return c[Symbol(ele)]
    end

    function Base.keys(c::Circuit)
        return c.keys
    end
