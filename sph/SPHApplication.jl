###############################################################################
## Application and basic Functions
###############################################################################
module SPHApplication
    # Import Solver
    using SPHSolver
    # Import Equation
    using SPHEquations
    # Import ParticlePool
    using SPHParticles

################################################################################
    abstract type AbstractApplication end
    # Basic Properties and general functions
    mutable struct Case <: AbstractApplication
        appsname::String
        foutdir::String
        solver::SPHSolver.Solver
        equations::SPHEquations.Group
        particles::SPHParticles.ParticlePoolUnion
        Case(aname::String="DefaultApplication",
            fdir::String="./DefaultApplication/") = new(aname,fdir)
    end

    # Casewise functions need to be specified in each case.
    function create_solver(p::SPHApplication.AbstractApplication) end
    function create_equations(p::SPHApplication.AbstractApplication) end
    function create_particles(p::SPHApplication.AbstractApplication) end

    function initialize(p::SPHApplication.AbstractApplication)
        start_time::Float64 = time()
        #
        p.particles = create_particles(p)
        p.equations = create_equations(p)
        p.solver = create_solver(p)
        SPHSolver.set_target(p.solver, p.equations, p.particles)
        SPHSolver.set_output_dir(p.solver, string("$(p.foutdir)/$(p.appsname)"))
        #
        end_time::Float64 = time()
        setup_duration::Float64 = end_time - start_time
        #write_info(info_filename, completed=False, cpu_time=0)
        display(p)
        println("Setup took: ","$setup_duration"," secs")
    end

    function apprun(p::SPHApplication.AbstractApplication)
        start_time::Float64 = time()
        SPHSolver.solve(p.solver)
        end_time::Float64 = time()
        run_duration::Float64 = end_time - start_time
        println("Run took: ","$(run_duration)"," secs")
        #write_info(info_filename, completed=True, cpu_time=run_duration)
    end

################################################################################
    function set_output_dir(p::SPHApplication.AbstractApplication, foutdir)
        p.foutdir = foutdir
    end

    function set_apps_name(p::SPHApplication.AbstractApplication, aname)
        p.appsname = aname
    end

    function Base.display(p::SPHApplication.AbstractApplication)
        println("***************************************************************")
        println("Case: $(p.appsname)        Output Directory: $(p.foutdir)")
        display(p.particles)
        display(p.equations)
        println("***************************************************************")
    end
################################################################################
end
