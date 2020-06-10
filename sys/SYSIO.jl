module SYSIO
    using JSON2
    using SYSEntity

##################################################################################
# Input Parts
##################################################################################
    function read_input(path::String)
        open(path, "r") do io
            return JSON2.read(join(readlines(io)))
        end
    end

##################################################################################
# Output Parts
##################################################################################
    mutable struct Output
        dir::String
        fd::IOStream
        stage::Int64
        stages::Vector{Tuple{Float64, Float64, Float64}}
        paras::Dict{Symbol,Vector{Any}}

        #
        hydro_nodes::Vector{SYSEntity.Node}
        thermal_nodes::Vector{SYSEntity.Node}

        Output(input::NamedTuple, dir::String) = init(new(), input, dir)
    end

    function init(opt::Output, input::NamedTuple, dir::String)
        if isdir(dir)
            rm(dir, recursive=true)
        end
        mkpath(dir)
        opt.dir = dir
        #
        opt.stage = 1.0
        opt.stages = Vector{Tuple{Float64, Float64, Float64}}()
        stage_num = length(input[:OUTPUT_TIME_SETTING][:TOTAL_TIME])-1
        for ii in 1:stage_num
            push!(opt.stages, (input[:OUTPUT_TIME_SETTING][:TOTAL_TIME][ii],
                           input[:OUTPUT_TIME_SETTING][:TOTAL_TIME][ii+1],
                           input[:OUTPUT_TIME_SETTING][:TIME_STEP][ii]))
        end
        #
        opt.paras = Dict{Symbol,Vector{Any}}()
        for key in keys(input[:OUTPUT_PARAMETER])
            push!(opt.paras, (key => input[:OUTPUT_PARAMETER][key]))
        end
        #
        return opt
    end

    function set_output(opt::Output, hydro::Vector{SYSEntity.Node}, thermal::Vector{SYSEntity.Node})
        opt.hydro_nodes = hydro
        opt.thermal_nodes = thermal
        output_header(opt)
    end

    function output_header(opt::Output)
        for key in [opt.paras[:PARAMETER];:pv]
            open(joinpath(opt.dir,String("hydro_output_$(String(key)).csv")), "a+") do io
                write(io, "TIME, ")
                for cv in opt.hydro_nodes
                    write(io, "$(cv.id), ")
                end
                write(io, "\n")
            end
        end
        open(joinpath(opt.dir,"thermal_output.csv"), "a+") do io
            write(io, "TIME, ")
            for grid in opt.thermal_nodes
                write(io, "$(grid.id), ")
            end
            write(io, "\n")
        end
    end

    function output(opt::Output, tc::Float64)
        set_stage(opt, tc)
        (ts, tf, dt) = opt.stages[opt.stage]
        err = (tc-ts)/dt%1.0
        if err < 1e-4 || err > 0.9999
            print("     output :         ")
            for key in opt.paras[:PARAMETER]
                open(joinpath(opt.dir, String("hydro_output_$(String(key)).csv")), "a+") do io
                    write(io, "$tc, ")
                    for cv in opt.hydro_nodes
                        write(io, "$(cv[key]), ")
                    end
                    write(io, "\n")
                end
            end
            open(joinpath(opt.dir,"hydro_output_pv.csv"), "a+") do io
                write(io, "$tc, ")
                for cv in opt.hydro_nodes
                    pv = cv[:p] + 0.5*cv[:rho]*cv[:v]^2.0
                    write(io, "$pv, ")
                end
                write(io, "\n")
            end
            open(joinpath(opt.dir,"thermal_output.csv"), "a+") do io
                write(io, "$tc, ")
                for grid in opt.thermal_nodes
                    write(io, "$(grid[:t]), ")
                end
                write(io, "\n")
            end
        end
    end

    function set_stage(opt::Output, tc::Float64)
        while tc > opt.stages[opt.stage][2]
            opt.stage += 1
        end
    end
end
