###############################################################################
## Application and basic Functions
###############################################################################
module SYSApplication
    using SYSIO
    using SYSEntity
    using SYSSolver
    using Controller

################################################################################
    abstract type AbstractApplication end
    # Basic Properties and general functions
    mutable struct Case <: AbstractApplication
        appsname::String
        inputpath::String
        foutdir::String
        #
        plant::SYSEntity.Plant
        solver::SYSSolver.Solver
        #
        tc::Controller.TimeSensor
        tcal::Vector{Tuple{Float64, Float64, Float64}}
        cons::Controller.ControlSystem
        #
        opt::SYSIO.Output
        Case(aname::String="DefaultSysApplication", input::String="Default.json",
            fdir::String="./DefaultSysApplication/") = new(aname, input, fdir)
    end

    function initialize(p::SYSApplication.AbstractApplication)
        start_time::Float64 = time()
        #----------------------------------------------------------------------------------
        input = SYSIO.read_input(p.inputpath)
        general_setting = input[:GENERAL_SETTING]
        # calculate time setting
        p.tc = Controller.TimeSensor(Float64.(general_setting[:CURRENT_TIME]))
        p.tcal = Vector{Tuple{Float64, Float64, Float64}}()
        stage_num = length(general_setting[:CAL_TIME_SETTING][:TOTAL_TIME])-1
        for ii in 1:stage_num
            push!(p.tcal, (general_setting[:CAL_TIME_SETTING][:TOTAL_TIME][ii],
                           general_setting[:CAL_TIME_SETTING][:TOTAL_TIME][ii+1],
                           general_setting[:CAL_TIME_SETTING][:TIME_STEP][ii]))
        end
        #----------------------------------------------------------------------------------
        # Plant Setup
        p.plant = SYSEntity.Plant(input[:PLANT_CONFIGURATION])
        #----------------------------------------------------------------------------------
        # Solver Setup
        p.solver = SYSSolver.Solver(p.plant, general_setting)
        #----------------------------------------------------------------------------------
        # Control System Setup
        p.cons = Controller.ControlSystem(input)
        Controller.add(p.cons, :TIME, p.tc)
        Controller.binding(p.cons, p.plant)
        #----------------------------------------------------------------------------------
        # Output Setup
        p.opt = SYSIO.Output(general_setting, joinpath(p.foutdir, p.appsname))
        SYSIO.set_output(p.opt, SYSEntity.get_all_cvs(p.plant), SYSEntity.get_all_grids(p.plant))
        #----------------------------------------------------------------------------------
        end_time::Float64 = time()
        setup_duration::Float64 = end_time - start_time
        display_info(p)
        println("Setup took: ","$setup_duration"," secs")
    end

    function apprun(p::SYSApplication.AbstractApplication)
        start_time::Float64 = time()
        for (ts, tf, dt) in p.tcal
            SYSSolver.set_dt(p.solver, dt)
            while p.tc.time <= tf
                println("t=$(p.tc.time)")
                Controller.action(p.cons)
                SYSSolver.solve(p.solver)
                @time SYSIO.output(p.opt, p.tc.time)
                p.tc.time += dt
            end
        end
        end_time::Float64 = time()
        run_duration::Float64 = end_time - start_time
        println("Run took: ","$(run_duration)"," secs")
        #write_info(info_filename, completed=True, cpu_time=run_duration)
    end

################################################################################
    function set_input_path(p::SYSApplication.AbstractApplication, input::String)
        p.input_file = input
    end

    function set_output_dir(p::SYSApplication.AbstractApplication, foutdir::String)
        p.foutdir = foutdir
    end

    function set_apps_name(p::SYSApplication.AbstractApplication, aname::String)
        p.appsname = aname
    end

    function display_info(p::SYSApplication.AbstractApplication)
        println("***********************************************************************************")
        println("Case: $(p.appsname)        Input File: $(p.inputpath)        Output Directory: $(p.foutdir)")
        SYSEntity.display_info(p.plant)
#        SYSSolver.display_info(p.solver)
        println("***********************************************************************************")
    end

    function output(p::SYSApplication.AbstractApplication)
        if mod(iteration_times, p.solver.pfreq)==0 || p.solver.tc in p.solver.output_at_times
            dirpath = string("$(p.foutdir)/$(p.appsname)/t=$(p.solver.tc)s")
        end
    end
################################################################################
end
