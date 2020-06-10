    ###############################################################################
    # Hydro nodes
    ###############################################################################
    mutable struct CV <: Node
        id::Int64
        ele::Element
        mat::MaterialLib.Material
        in_ele_id::Int64
        #
        data::Dict{Symbol,Float64}
        links::Vector{Tuple{Face, Node, Int64}}
        #
        CV() = new(0)
        CV(id::Int64, ele::Element, in_ele_id::Int64, fluid::MaterialLib.Material) =
        new(id, ele, fluid, in_ele_id,
                    Dict{Symbol,Float64}(),
                    Vector{Tuple{Face, Node, Int64}}())
    end

    mutable struct TDV <: Node
        id::Int64
        ele::Element
        mat::MaterialLib.Material
        in_ele_id::Int64
        #
        data::Dict{Symbol,Float64}
        links::Vector{Tuple{Face, Node, Int64}}
        TDV(id::Int64, ele::Element, fluid::MaterialLib.Material) =
        new(id, ele, fluid, 1,
                    Dict{Symbol,Float64}(),
                    Vector{Tuple{Face, Node, Int64}}())
    end

    #------------------------------------------------------------------------------
    mutable struct Junction{CV_PREV<:Node,CV_NEXT<:Node} <: Face
        id::Int64
        ele::Element
        #
        data::Dict{Symbol,Float64}
        prev::CV_PREV
        next::CV_NEXT
        Junction(id::Int64, ele::Element, prev, next) =
        new{typeof(prev), typeof(next)}(id, ele, Dict{Symbol,Float64}(), prev, next)
    end

    mutable struct TDJ{CV_PREV<:Node,CV_NEXT<:Node} <: Face
        id::Int64
        ele::Element
        #
        data::Dict{Symbol,Float64}
        prev::CV_PREV
        next::CV_NEXT
        TDJ(id::Int64, ele::Element, prev, next) =
        new{typeof(prev), typeof(next)}(id, ele, Dict{Symbol,Float64}(), prev, next)
    end

    #------------------------------------------------------------------------------
    function init(cv::Node, fluid::MaterialLib.Material)
        temperature = cv.ele[:TEMPERATURE][1]
        if cv.ele[:NODE_NUM] >1
            temperature += (cv.ele[:TEMPERATURE][end]-cv.ele[:TEMPERATURE][1])*(cv.in_ele_id-1)/(cv.ele[:NODE_NUM]-1)
        end
        push!(cv.data, (:t => temperature))
        push!(cv.data, (:p => cv.ele[:PRESSURE][cv.in_ele_id]))
        push!(cv.data, (:p_n => cv.ele[:PRESSURE][cv.in_ele_id]))
        push!(cv.data, (:rho => MaterialLib.ro(fluid, cv[:t], cv[:p])))
        push!(cv.data, (:rho_n => MaterialLib.ro(fluid, cv[:t], cv[:p])))
        push!(cv.data, (:enthalpy => MaterialLib.enth(fluid, cv[:t], cv[:p])))
        push!(cv.data, (:enthalpy_n => MaterialLib.enth(fluid, cv[:t], cv[:p])))
        velocity::Float64 = cv.ele[:FLOWRATE][cv.in_ele_id] / cv[:rho] / cv.ele[:AREA]
        push!(cv.data, (:v => velocity))
        push!(cv.data, (:flow_resistance => 0.0))
        push!(cv.data, (:h_coeff => 0.0))
        push!(cv.data, (:lamda => MaterialLib.lamda(fluid, cv[:t], cv[:p])))
        push!(cv.data, (:cp => MaterialLib.c_p(fluid, cv[:t], cv[:p])))
        push!(cv.data, (:q_bak => cv.ele[:HEATSOURCE][cv.in_ele_id]))
        push!(cv.data, (:q => cv.ele[:HEATSOURCE][cv.in_ele_id]))
    end

    function init(junc::Face, fluid::MaterialLib.Material)
        push!(junc.data, (:rho => junc.prev[:rho]))
        push!(junc.data, (:enthalpy => junc.prev[:enthalpy]))
        push!(junc.data, (:p => junc.prev[:p]))
        push!(junc.data, (:v => junc.prev[:v]))
        push!(junc.data, (:cv_n => junc[:v]))
        push!(junc.data, (:head_source_a => 0.0))
        push!(junc.data, (:area => min(junc.prev.ele[:AREA], junc.next.ele[:AREA])))
    end

    function link(junc::Face)
        # add links
        push!(junc.prev.links, (junc, junc.next, -1))
        push!(junc.next.links, (junc, junc.prev, +1))
    end

    ################################################################################
    function set_temperature(cv::Node, temperature::Float64)
        #println("BEFORE   $(cv[:t])   $(cv[:rho])   $(cv[:enthalpy])")
        cv[:t] = temperature
        cv[:rho] = MaterialLib.ro(cv.mat, cv[:t], cv[:p])
        cv[:enthalpy] = MaterialLib.enth(cv.mat, cv[:t], cv[:p])
        cv.links[1][1][:rho] = cv[:rho]
        cv.links[1][1][:enthalpy] = cv[:enthalpy]
        #println("AFTER    $(cv[:t])   $(cv[:rho])   $(cv[:enthalpy])")
    end

    function set_temperature(cv::TDV, temperature::Nothing)
        nothing
    end

    function set_flowrate(junc::TDJ, flowrate::Float64)
        junc[:v] = flowrate / junc[:area] / junc[:rho]
        junc[:cv_n] = junc[:v]
        junc.prev[:v] = junc[:v]
    end
