    ###############################################################################
    abstract type AbstractSensor end
    function signal(s::AbstractSensor) nothing end

    ###############################################################################
    mutable struct TimeSensor <: AbstractSensor
        time::Float64
    end

    function signal(s::TimeSensor)
        return s.time
    end

    function update(s::TimeSensor, dt::Float64)
        s.time += dt
    end

    ###############################################################################
    mutable struct FlowrateSensor <: AbstractSensor
        cv::SYSEntity.CV
    end

    function signal(s::FlowrateSensor)
        return s.cv[:v] * s.cv[:rho] * s.cv.ele[:AREA]
    end

    ###############################################################################
    mutable struct TemperatureSensor{Type<:SYSEntity.Node} <: AbstractSensor
        node::Type
    end

    function signal(s::TemperatureSensor)
        return s.node[:t]
    end

    ###############################################################################
    mutable struct PointSensor{node <: SYSEntity.Node} <: AbstractSensor
        variable::Symbol
        thresholds::Vector{Float64}
        logic::Expr
        watcher::node
    end

    function signal(s::PointSensor)
        #ex = Expr(:call, logic, $watcher[variable], $threshold)
        #return eval(ex)
        nothing
    end

    ##########################################################################
    mutable struct ElementSensor{node <: SYSEntity.Node} <: AbstractSensor
        variable::Symbol
        thresholds::Vector{Float64}
        logic::Expr
        watcher::Vector{node}
    end
