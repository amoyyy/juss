    mutable struct Hydro1DMesh <: Mesh
        circuit::Circuit
        mat::MaterialLib.Material
        cvs::Vector{Node}
        juncs::Vector{Face}
        #
        _start_cv::Dict{Symbol,Int64}
        _end_cv::Dict{Symbol,Int64}

        Hydro1DMesh(circuit::Circuit) =
        init(new(circuit, MaterialLib.Sodium(),
                          Vector{Node}(), Vector{Face}(),
                          Dict{Symbol,Int64}(), Dict{Symbol,Int64}()))
    end

    function init(s::Hydro1DMesh)
        s.mat = MaterialLib.getMat(s.circuit.fluid)
        # add node
        for key in keys(s.circuit)
            ele = s.circuit[key]
            start_cv::Int64 = length(s.cvs) + 1
            push!(s._start_cv, (key => start_cv))
            if is_type(ele, [:TDV, :TDJ])
                add_node(s, TDV(start_cv, ele, s.mat))
            else
                foreach(ii->add_node(s, CV(start_cv+ii-1, ele, ii, s.mat)), 1:ele[:NODE_NUM])
                foreach(ii->add_node(s, Junction(length(s.juncs)+1, ele, s.cvs[start_cv+ii-1], s.cvs[start_cv+ii])), 1:ele[:NODE_NUM]-1)
            end
            push!(s._end_cv, (key => length(s.cvs)))
        end
        for key in keys(s.circuit)
            ele = s.circuit[key]
            if is_type(ele, [:TDJ])
                add_node(s, TDJ(length(s.juncs)+1, ele, s.cvs[s._end_cv[key]], s.cvs[s._start_cv[ele.next[1].name]]))
            elseif is_type(ele, [:VALVE])
                nothing
            elseif ele[:NODE_NUM]>1
                foreach(next->add_node(s, Junction(length(s.juncs)+1, ele, s.cvs[s._end_cv[key]], s.cvs[s._start_cv[next.name]])),ele.next)
            end
        end
        # init and link node
        foreach(cv->init(cv, s.mat), s.cvs)
        foreach(junc->init(junc, s.mat), s.juncs)
        foreach(link, s.juncs)
        #
        return s
    end

    function add_node(s::Hydro1DMesh, face::Face)
        push!(s.juncs, face)
    end

    function add_node(s::Hydro1DMesh, node::Node)
        push!(s.cvs, node)
    end

    function Base.getindex(s::Hydro1DMesh, ele::Symbol)
        return s.cvs[s._start_cv[ele]:s._end_cv[ele]]
    end

    function Base.getindex(s::Hydro1DMesh)
        return s.cvs[1:end]
    end
