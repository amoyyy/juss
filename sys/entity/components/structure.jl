    ##------------------------------------------------------------------------------
    ## Structure
    ##------------------------------------------------------------------------------
    mutable struct Structure <: Component
        name::Symbol
        keys::Dict{Symbol, DataType}
        data::Dict{DataType, Dict{Symbol, Any}}
        #

        Structure(name::Symbol) = new(name, Dict{Symbol, DataType}(), Dict{DataType, Dict{Symbol, Any}}())
        Structure(name::Symbol, input::NamedTuple) = init(Structure(name), input)
    end

    function init(t::Structure, input::NamedTuple)
        #####################################################################
        INT_KEY = [:NODE_NUM_AXIAL, :NODE_NUM_RADIAL, :INNER_BOUNDARY_TYPE, :OUTER_BOUNDARY_TYPE]
        SYMBOL_KEY = [:MATERIAL_MEDIA, :INNER_BOUNDARY_LINK, :OUTER_BOUNDARY_LINK]
        FLOAT_KEY = [:LENGTH]
        VECOTR_KEY_1 = [:INNER_BOUNDARY_VALUE, :OUTER_BOUNDARY_VALUE]
        VECOTR_KEY_2 = [:TEMPERATURE]
        #####################################################################
        foreach(key->add_value(t, key, Int64(input[key])), INT_KEY)
        foreach(key->add_value(t, key, Symbol(input[key])), SYMBOL_KEY)
        foreach(key->add_value(t, key, Float64(input[key])), FLOAT_KEY)
        for key in VECOTR_KEY_1
            data_vector::Vector{Float64} = Array(input[key])
            last_value::Float64 = data_vector[length(data_vector)]
            foreach(ii->push!(data_vector, last_value), 1:t[:NODE_NUM_AXIAL]-length(data_vector))
            add_value(t, key, data_vector)
        end
        for key in VECOTR_KEY_2
            data_vector::Vector{Float64} = Array(input[key])
            last_value::Float64 = data_vector[length(data_vector)]
            foreach(ii->push!(data_vector, last_value), 1:t[:NODE_NUM_AXIAL]*t[:NODE_NUM_RADIAL]-length(data_vector))
            add_value(t, key, data_vector)
        end
        add_value(t, :INNER_CORDIN, input[:COORDINATE_LAYER][1])
        add_value(t, :OUTER_CORDIN, input[:COORDINATE_LAYER][2])
        add_value(t, :DELTA_X, t[:LENGTH]/t[:NODE_NUM_AXIAL])
        add_value(t, :DELTA_R, (t[:OUTER_CORDIN]-t[:INNER_CORDIN])/t[:NODE_NUM_RADIAL])
        return t
    end

    function Base.size(s::Structure)
        return tuple(s[:NODE_NUM_AXIAL], s[:NODE_NUM_RADIAL])
    end
