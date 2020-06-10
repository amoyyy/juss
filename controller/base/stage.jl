    mutable struct Stage
        t_start::Float64
        t_end::Float64
        dt::Float64
        Stage(t_start::Float64, t_end::Float64, dt::Float64=0.0) = new(t_start, t_end, dt)
    end

    mutable struct TimeStages
        stage::Int
        time_stages::Vector{Stage}
        TimeStages() = new(1, Vector{Stage}())
        TimeStages(time_input::NamedTuple) = init(TimeStages(), time_input)
    end

    function init(con::TimeStages, time_input::NamedTuple)
        tfs = Float64.(time_input[:TOTAL_TIME])
        dts = Float64.(time_input[:TIME_STEP])
        #@assert length(tfs) == length(dts) "ERROR: Time Setting Error, $time_input"
        for idx in 1:length(tfs)-1
            push!(con.time_stages, Stage(tfs[idx], tfs[idx+1], dts[idx]))
        end
        return con
    end

    stage_num = (s::TimeStages) -> length(s.time_stages)
    ts = (s::TimeStages) -> s.time_stages[s.stage].t_start
    tf = (s::TimeStages) -> s.time_stages[s.stage].t_end
    dt = (s::TimeStages) -> s.time_stages[s.stage].dt

    function set_stage(con::TimeStages, tc::Float64)
        while tc > tf(con)
            con.stage += 1
            if(con.stage>length(con.time_stages))
                println("ERROR:  stage error!")
            end
        end
    end
