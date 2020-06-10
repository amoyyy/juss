    ##------------------------------------------------------------------------------
    ## Element
    ##------------------------------------------------------------------------------
    mutable struct Element <: Component
        name::Symbol
        keys::Dict{Symbol, DataType}
        data::Dict{DataType, Dict{Symbol, Any}}
        prev::Vector{Element}
        next::Vector{Element}
        Element(name::Symbol) = new(name,
                                    Dict{Symbol, DataType}(),
                                    Dict{DataType, Dict{Symbol, Any}}(),
                                    Vector{Element}(),
                                    Vector{Element}())
        Element(name::Symbol, input::NamedTuple) = init(Element(name), input)
    end

    function init(e::Element, input::NamedTuple)
        #####################################################################
        INT_KEY = [:MULTI_COPY, :NODE_NUM]
        FLOAT_KEY = [:LENGTH, :AREA, :WETTED_PERIMETER, :THETA, :PHI]
        VECOTR_KEY = [:FLOWRATE, :TEMPERATURE, :PRESSURE, :HEATSOURCE, :ENTHALPY_RESTART]
        #####################################################################
        foreach(key->add_value(e, key, Int64(input[key])), INT_KEY)
        foreach(key->add_value(e, key, Float64(input[key])), FLOAT_KEY)
        for key in VECOTR_KEY
            data_vector::Vector{Float64} = Array(input[key])
            last_value::Float64 = data_vector[length(data_vector)]
            foreach(ii->push!(data_vector, last_value), 1:e[:NODE_NUM]-length(data_vector))
            add_value(e, key, data_vector)
        end
        add_value(e, :DELTA_X, e[:LENGTH]/e[:NODE_NUM])
        add_value(e, :HYDRO_DIAMETER, 4.0*e[:AREA]/e[:WETTED_PERIMETER])
        return e
    end

    function add_next(p::Element, next::Element)
        push!(p.next, next)
        push!(next.prev, p)
    end

    #---------------------------------------------------------------------------------
    function is_type(e::Element, type::Symbol)
        # judge whether an Element is a specific type according to the data...
        conditon_1 = occursin(string(type), string(e.name))   # name condition .. not precise, promoting in future
        #condition_2 = true ....
        if conditon_1
            return true
        end
        return false
    end

    function is_type(e::Element, types::Vector{Symbol})
        for type in types
            if is_type(e, type)
                return true
            end
        end
        return false
    end
