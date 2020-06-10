###############################################################################
## Controller
###############################################################################
module Controller
    using SYSEntity
    using SYSSolver
#    using SPHParticles
#    using NNP

    include("base/funcs.jl")
    include("input/curves.jl")
    include("input/tables.jl")

    include("sensor/sensor.jl")
    include("controller/controller.jl")

    mutable struct ControlSystem
        sensors::Dict{Symbol, AbstractSensor}
        controllers::Dict{Symbol, AbstractController}
        paras::Dict{Symbol, Any}

        ControlSystem(input::NamedTuple) =
        new(Dict{Symbol, AbstractSensor}(), Dict{Symbol, AbstractController}(),
            Dict{Symbol, Any}())
    end

    function action(c::ControlSystem)
        foreach(action, values(c.controllers))
    end

    function add(c::ControlSystem, key::Symbol, s::AbstractSensor)
        push!(c.sensors, (key => s))
    end

    function add(c::ControlSystem, key::Symbol, con::AbstractController)
        push!(c.controllers, (key => con))
    end

    function sensor(c::ControlSystem, key::Symbol)
        if key in keys(c.sensors)
            return c.sensors[key]
        end
        nothing
    end

    function controller(c::ControlSystem, key::Symbol)
        if key in keys(c.controllers)
            return c.controllers[key]
        end
        nothing
    end

    ##########################################################################
    function binding(c::ControlSystem, s::SYSSolver.Solver)
    end

    function binding(c::ControlSystem, p::SYSEntity.Plant)
        tdjv01 = SYSEntity.get_cvs(p, :TDJ01)[1]
        tdjv02 = SYSEntity.get_cvs(p, :TDJ02)[1]
        # ------------------------------------------------
        # case1 CFER IHX
        #add(c, :TDJ01_FLOWRATE, InletController(sensor(c,:TIME), tdjv01, flowrate_curve_1, temperature_constant))
        #add(c, :TDJ02_FLOWRATE, InletController(sensor(c,:TIME), tdjv02, flowrate_curve_2, temperature_constant))
        # ------------------------------------------------
        # case2 CFR600 IHX
        add(c, :TDJ01_FLOWRATE, InletController(sensor(c,:TIME), tdjv01, flowrate_table_1, temperature_table_1))
        add(c, :TDJ02_FLOWRATE, InletController(sensor(c,:TIME), tdjv02, flowrate_table_2, temperature_table_2))
        # ------------------------------------------------
        #add(c, :TDJ01_FLOWRATE, InletController(sensor(c,:TIME), tdjv01, flowrate_sind, temperature_constant))
    end

end
