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
    L = 1.0
    nx = 50
    dx = L/nx
    extention = 5*dx
    volume = dx*dx

    gy = 0.0
    vmax = 1.0
    c0 = 10 * vmax
    re = 100
    nu = vmax*L/re

    rho0 = 1.0
    p0 = c0*c0*rho0
    gamma = 1.0

    VARS = [:x, :y, :z, :h, :p, :rho, :m, :V, :T, :H, :aT, :vx, :vy, :vz, :ax, :ay, :az, :arho, :cs]
    VARS_FLUID_ADDITION = [:vx0, :vy0, :vz0, :ax0, :ay0, :az0]
    VARS_BOUNDARY_ADDITION = [:nx, :ny, :nz, :tx, :ty, :tz]

    VARSALL = [:x, :y, :z, :h, :p, :rho, :m, :V, :T, :H, :aT, :vx, :vy, :vz, :ax, :ay, :az, :arho,
               :rho0, :V0, :x0, :y0, :z0, :vx0, :vy0, :vz0, :ax0, :ay0, :az0, :vxf, :vyf, :vzf, :vxg, :vyg, :vzg, :wij]
###############################################################################
    # DEBUG
    function SPHApplication.create_solver(p::SPHApplication.Case)
        dt_cfl = 0.25 * h0/( c0 + vmax )
        dt_viscous = 0.125 * h0^(2)/nu
        dt_force = 1.0
        dt = min(dt_cfl, dt_viscous, dt_force)
        tf = 10.0
        solver = SPHSolver.Solver(2,
                        #SPHIntegrators.PECIntegrator(:Fluid, SPHIntegrators.TransportVelocityStep()),
                        SPHIntegrators.PECIntegrator(:Fluid, SPHIntegrators.WCSPHStep()),
                        #QuinticSpline(2),
                        SPHKernels.QuinticSpline(2),
                        tf, 0.0, dt, hdx, 100, [0.0, tf])
        return solver
    end

    function SPHApplication.create_equations(p::SPHApplication.Case)
        equations = SPHEquations.Group()
"""        SPHEquations.add(equations, 1, SPHEquations.NumberDensity(:Fluid, :ALL))
        SPHEquations.add(equations, 2, SPHEquations.TaitEOS(:Fluid, :NONE, rho0, c0, 0.0, gamma))
        SPHEquations.add(equations, 2, SPHEquations.SetWallVelocity(:Solid, :Fluid))
        SPHEquations.add(equations, 3, SPHEquations.SolidWallPressureBC(:Solid, :Fluid, rho0, p0, 0.0, 0.0, 0.0))
        SPHEquations.add(equations, 4, SPHEquations.MomentumEquationPressureGradient(:Fluid, :ALL, p0, 0.0, 0.0, 0.0))
        #SPHEquations.add(equations, 4, SPHEquations.ArtificialViscosity(:Fluid, :ALL, 0.01, c0))
        SPHEquations.add(equations, 4, SPHEquations.Viscosity(:Fluid, :Fluid, nu))
        SPHEquations.add(equations, 4, SPHEquations.WallNoSlipBC(:Fluid, :Solid, nu))
        SPHEquations.add(equations, 4, SPHEquations.ArtificialStress(:Fluid, :Fluid))"""
        SPHEquations.add(equations, 2, SPHEquations.VolumeFromMassDensity(:Fluid))
        SPHEquations.add(equations, 2, SPHEquations.TaitEOS(:Fluid, :NONE, rho0, c0, 0.0, gamma))
        SPHEquations.add(equations, 3, SPHEquations.SetWallVelocity(:Solid, :Fluid))
        SPHEquations.add(equations, 3, SPHEquations.SolidWallPressureBC(:Solid, :Fluid, rho0, p0, 0.0, 0.0, 0.0))
        SPHEquations.add(equations, 4, SPHEquations.ContinuityEquation(:Fluid, :ALL))
        SPHEquations.add(equations, 4, SPHEquations.MomentumEquationPressureGradient(:Fluid, :ALL, p0, 0.0, 0.0, 0.0))
        SPHEquations.add(equations, 4, SPHEquations.ArtificialViscosity(:Fluid, :ALL, 1.0, c0))
        #SPHEquations.add(equations, 4, SPHEquations.Viscosity(:Fluid, :Fluid, nu))
        SPHEquations.add(equations, 4, SPHEquations.WallNoSlipBC(:Fluid, :Solid, nu))
        #SPHEquations.add(equations, 4, SPHEquations.XSPHCorrection(:Fluid, :Fluid, 0.1))
        return equations
    end

    function SPHApplication.create_particles(p::SPHApplication.Case)
        p = SPHParticles.ParticlePoolUnion()
        p.comments="
        #      ______________________
        #      _________s ___________
        #   | |                      | |
        #   | |                      | |
        #   |s|          f           |s|
        #   | |                      | |
        #   | |                      | |
        #   | |______________________| |
        #      _________s____________
        "
        # Add Pools
        RANGE = [-0.11,1.11]
        NDIMS = [62,62]
        x0 = GeoLib.regular(RANGE, RANGE, NDIMS)[1]
        y0 = GeoLib.regular(RANGE, RANGE, NDIMS)[2]
        fluid_idx = (x0.>0.0).&(x0.<1.0).&(y0.>0.0).&(y0.<1.0)
        solid_idx = .!fluid_idx
        x_fluid = x0[fluid_idx]
        y_fluid = y0[fluid_idx]
        x_solid = x0[solid_idx]
        y_solid = y0[solid_idx]
        TAG = [:Fluid, :Solid]
        NVP = [length(x_fluid), length(x_solid)]

        for ii in 1:length(TAG)
            pool = SPHParticles.ParticlePool(ii, TAG[ii], NVP[ii])
            SPHParticles.set_var(pool, :m, rho0*volume)
            SPHParticles.set_var(pool, :h, hdx*dx)
            if TAG[ii] == :Solid
                v = zeros(NVP[ii])
                v[y_solid.>1.0] .= 1.0
                SPHParticles.set_var(pool, :vx, v)
                SPHParticles.set_var(pool, :x, x_solid)
                SPHParticles.set_var(pool, :y, y_solid)
                SPHParticles.set_var(pool, :rho, rho0)
                SPHParticles.set_var(pool, :V, 1.0/volume)
            elseif TAG[ii] == :Fluid
                rho_modify = ((rand(NVP[ii]).-0.5)/500 .+1.0)*rho0
                SPHParticles.set_var(pool, :rho, rho_modify)
                SPHParticles.set_var(pool, :V, rho_modify./volume./rho0)
                SPHParticles.set_var(pool, :x, x_fluid)
                SPHParticles.set_var(pool, :y, y_fluid)
            end
            SPHParticles.add(p, TAG[ii], pool)
        end

        # Set Boundarys
        SPHParticles.link(p[:Solid,1], p[:Fluid,1] ,GeoLib.Boundary(GeoLib.NoBoundary(), ==))
        return p
    end

    hdx = 1.0
    h0 = hdx * dx
    for hdx_tmp in [1.0, 1.6, 2.0, 3.1]
        hdx = hdx_tmp
        h0 = hdx * dx
        case = SPHApplication.Case("hdx=$hdx", "./cases/output/Cavity3")
        SPHApplication.initialize(case)
        SPHApplication.apprun(case)
    end
