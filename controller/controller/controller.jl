    ###############################################################################
    abstract type AbstractInput end
    mutable struct Curve <: AbstractInput
        variable::Vector{Symbol}
        curve::Expr

        Curve() = new(Vector{Symbol}(), Expr(:NONE))
        Curve(variables::Vector{Symbol}, func::Expr) = new(variables, func)
    end

    mutable struct Table <: AbstractInput
        variable::Symbol
        table::Tuple{Pair{Float64, Float64}, Float64}
    end
    ###############################################################################
    abstract type AbstractController end
    function binding(con::AbstractController, sensor::AbstractSensor) end
    function action(con::AbstractController) end

    ###############################################################################
    mutable struct FlowrateController <: AbstractController
        sensor::AbstractSensor
        exprs::Function
        object::SYSEntity.Face

        FlowrateController(s::AbstractSensor, ex::Function, object::SYSEntity.Face) =
        new(s, ex, object)
        #CurveController(time_input::NamedTuple, paras_input::NamedTuple) =
        #init(new(), time_input, paras_input)
    end

    function action(c::FlowrateController)
        factor = c.exprs(signal(c.sensor))
        #println("$(signal(c.sensor))   $factor")
        flowrate = c.object.ele[:FLOWRATE][1]*factor
        SYSEntity.set_flowrate(c.object, flowrate)
    end

    ###############################################################################
    mutable struct TemperatureController <: AbstractController
        sensor::AbstractSensor
        exprs::Function
        object::SYSEntity.Node

        TemperatureController(s::AbstractSensor, ex::Function, object::SYSEntity.Node) =
        new(s, ex, object)
        #CurveController(time_input::NamedTuple, paras_input::NamedTuple) =
        #init(new(), time_input, paras_input)
    end

    function action(c::TemperatureController)
        temperature = c.exprs(signal(c.sensor))
        #println(" $(signal(c.sensor))    $temperature")
        SYSEntity.set_temperature(c.object, temperature)
    end

    ###############################################################################
    mutable struct InletController <: AbstractController
        sensor::AbstractSensor
        object1::SYSEntity.Node
        object2::SYSEntity.Face
        flowrate::Function
        temperature::Function

        InletController(s::AbstractSensor, object::SYSEntity.Node, ex1::Function, ex2::Function) =
        new(s, object, object.links[1][1], ex1, ex2)
    end

    function action(c::InletController)
        time = signal(c.sensor) - 10.0
        if time>=0
            flowrate = c.flowrate(time)*c.object1.ele[:FLOWRATE][1]
            temperature = c.temperature(time)
            #println(time,"   ",flowrate, "   ",temperature, "   ",c.object1[:t])
            SYSEntity.set_temperature(c.object1, temperature)
            SYSEntity.set_flowrate(c.object2, flowrate)
        end
    end
