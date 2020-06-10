include("../Dependency.jl")
# Common modules used
using MaterialLib
using GeoLib
# Specific modules that must be used
using SPHParticles
using SPHKernels
using SPHEquations
using SPHIntegrators
using SPHSolver
using SPHApplication

###############################################################################
## Several Conditions
###############################################################################
    L = 0.04
    nx = 20
    Umax = 1.0
    c0 = 10 * Umax
    rho0 = 860.0
    p0 = 100000.0
    dx = L/nx
    volume = (dx*dx)
    hdx = 1.6
    T0 = 500.0
    H = MaterialLib.enth(MaterialLib.Sodium(), T0, p0)
"""    re = 100
    h0 = hdx * dx
    nu = Umax*L/re
    dt_cfl = 0.25 * h0/( c0 + Umax )
    dt_viscous = 0.125 * h0^(2)/nu
    dt_force = 1.0
    dt = min(dt_cfl, dt_viscous, dt_force)"""
    dt = 1e-2
    tf = 5.0

    VARSALL = [:x, :y, :z, :h, :p, :rho, :m, :V, :T, :H,
               :aT, :aH, :vx, :vy, :vz, :ax, :ay, :az, :arho]
###############################################################################
    # DEBUG
    function SPHApplication.create_solver(p::SPHApplication.Case)
        solver = SPHSolver.Solver(2,
                        #PECIntegrator(:Fluid, TransportVelocityStep()),
                        SPHIntegrators.EulerIntegrator(:Fluid, SPHIntegrators.EulerStep()),
                        #QuinticSpline(2),
                        SPHKernels.Gaussian(2),
                        tf, 0.0, dt, hdx,
                        100, [0.0, tf])
        return solver
    end

    function SPHApplication.create_equations(p::SPHApplication.Case)
        equations = SPHEquations.Group()
        SPHEquations.add(equations, 1, SPHEquations.NumberDensity(:Fluid, :ALL))
        SPHEquations.add(equations, 2, SPHEquations.EnergyEquationH(:Fluid, :ALL, MaterialLib.Sodium(), MaterialLib.Sodium(), 0.00))
        return equations
    end

    function SPHApplication.create_particles(p::SPHApplication.Case)
        p = SPHParticles.ParticlePoolUnion()
        p.comments="
        #
        #      ______________________
        #   | |                      | |
        #   | |                      | |
        #   |s|          f1          |s|
        #   |1|                      |2|
        #   | |                      | |
        #   | |______________________| |
        "
        # Add Pools
        NUM_FLUID::Int = 400
        NUM_SOLID::Int = 100

        XRANGE_FLUID = [[0.0, 0.04]]
        YRANGE_FLUID = [[0.0, 0.04]]

        XRANGE_SOLID = [[-0.00,0.00],[0.04,0.04]]
        YRANGE_SOLID = [[0.0, 0.04],[0.0, 0.04]]
        NDIMS_SOLID = [[1,100], [1,100]]
        T_SOLID = [250.0, 250.0]
        H_SOLID = [MaterialLib.enth(MaterialLib.Sodium(), 250.0, p0), MaterialLib.enth(MaterialLib.Sodium(), 250.0, p0)]

        for ii in Vector(1:1)
            fluid = SPHParticles.ParticlePool(ii, :Fluid, NUM_FLUID)
            SPHParticles.set_var(fluid, :x, GeoLib.regular(XRANGE_FLUID[ii], YRANGE_FLUID[ii], [20,20])[1])
            SPHParticles.set_var(fluid, :y, GeoLib.regular(XRANGE_FLUID[ii], YRANGE_FLUID[ii], [20,20])[2])
            SPHParticles.set_var(fluid, :rho, rho0)
            SPHParticles.set_var(fluid, :V, 1.0/volume)
            SPHParticles.set_var(fluid, :m, rho0*volume)
            SPHParticles.set_var(fluid, :h, hdx*dx)
            SPHParticles.set_var(fluid, :p, p0)
            SPHParticles.set_var(fluid, :T, T0)
            SPHParticles.set_var(fluid, :H, H)
            SPHParticles.add(p, :Fluid, fluid)
        end
        for ii in Vector(1:2)
            solid = SPHParticles.ParticlePool(ii, :Solid, NUM_SOLID)
            SPHParticles.set_var(solid, :x, GeoLib.regular(XRANGE_SOLID[ii], YRANGE_SOLID[ii], NDIMS_SOLID[ii])[1])
            SPHParticles.set_var(solid, :y, GeoLib.regular(XRANGE_SOLID[ii], YRANGE_SOLID[ii], NDIMS_SOLID[ii])[2])
            SPHParticles.set_var(solid, :rho, rho0)
            SPHParticles.set_var(solid, :V, 1.0/volume)
            SPHParticles.set_var(solid, :m, rho0*volume)
            SPHParticles.set_var(solid, :h, hdx*dx)
            SPHParticles.set_var(solid, :p, p0)
            SPHParticles.set_var(solid, :T, T_SOLID[ii])
            SPHParticles.set_var(solid, :H, H_SOLID[ii])
            SPHParticles.add(p, :Solid, solid)
        end

        # Set Boundarys
        boundarys = GeoLib.Boundary[]
        push!(boundarys, GeoLib.Boundary(GeoLib.PlateX(0.00), <))
        push!(boundarys, GeoLib.Boundary(GeoLib.PlateX(0.04), <))
        SPHParticles.link([[p[:Solid,1]], [p[:Fluid,1]], [p[:Solid,2]]], boundarys)
        return p
    end


case = SPHApplication.Case("SPH-HeatConduction_H", "./cases/output")
SPHApplication.initialize(case)
SPHApplication.apprun(case)
