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
    H = 0.32
    P = 0.08
    L = 0.16
    dx = dy = dz = 0.01
    volume = dx*dy*dz

    gy = -9.8
    Vmax = sqrt(abs(gy) * H)
    c0 = 10 * Vmax

    rho0 = 860.0
    p0 = c0*c0*rho0

    hdx = 1.6
    h0 = hdx * dx
    re = 100.0
    nu = Vmax*P/re

    dt_cfl = 0.25 * h0/( c0 + Vmax )
    dt_viscous = 0.125 * h0^(2)/nu
    dt_force = 1.0
    dt = min(dt_cfl, dt_viscous, dt_force)*10
    tf = 1.0

    VARS = [:x, :y, :z, :h, :p, :rho, :m, :V, :T, :H, :aT, :vx, :vy, :vz, :ax, :ay, :az, :arho, :cs]
    VARS_FLUID_ADDITION = [:x_bak, :y_bak, :vx_bak, :vy_bak, :rho_bak]
    VARS_BOUNDARY_ADDITION = [:nx, :ny, :nz, :tx, :ty, :tz]
###############################################################################
    # DEBUG
    function SPHApplication.create_solver(p::SPHApplication.Case)
        solver = SPHSolver.Solver(3,
                        #PECIntegrator(:Fluid, TransportVelocityStep()),
                        SPHIntegrators.PECIntegrator(:Fluid, SPHIntegrators.WCSPHStep()),
                        #QuinticSpline(2),
                        SPHKernels.CubicSpline(3),
                        tf, 0.0, dt, hdx, 1, [0.0, tf])
        return solver
    end

    function SPHApplication.create_equations(p::SPHApplication.Case)
        equations = SPHEquations.Equation[]
        push!(equations, SPHEquations.NumberDensityEquation(:Fluid, :Fluid))
        push!(equations, SPHEquations.TaitEOS(:Fluid, 1e3, 10.0, 1e5, 7.0))
        push!(equations, SPHEquations.MomentumEquation(:Fluid, c0, 0.01, 0.0, 0.0, 0.0, 0.0))
        #push!(equations, SPHEquations.MonaghanBoundaryForce(0.01))
        push!(equations, SPHEquations.MonaghanKajtarBoundaryForce(:Fluid, :Solid, 0.02, 1.0, h0))
        return SPHEquations.Group(equations)
    end

    function SPHApplication.create_particles(p::SPHApplication.Case)
        p = SPHParticles.ParticlePoolUnion()
        p.comments="
        #          ______________________
        #        /|         P=0.8       /|
        #  L=1.6/ |                    / |
        #      /  |          4        /  |H=3.2
        #     /___|__________________/   |
        #     |   |__________________|___|
        #     | 1 /                  | 3 /
        #     |  /       2           |  /
        #     | /            5       | /
        #     |/_____________________|/
        #
        #            y|  /z
        #             | /
        #             |/----> x
        "
        # Add Pools
        NUM_POOL = 6
        TAG = [:Fluid, :Solid, :Solid, :Solid, :Solid, :Solid]
        NUM = [2560, 10240, 5120, 10240, 5120, 2560]
        CORDIN = [[[0.0, 0.08], [0.0, 0.20], [0.0, 0.16]],
                  [[-0.002, 0.0], [0.0, 0.32], [-0.002, 0.162]],
                  [[-0.002, 0.082], [0.0, 0.32], [-0.002, 0.0]],
                  [[0.08, 0.082], [0.0, 0.32], [-0.002, 0.162]],
                  [[-0.002, 0.082], [0.0, 0.32], [0.160, 0.162]],
                  [[-0.002, 0.082], [-0.002, 0.0], [-0.002, 0.162]]]
        NDIMS = [[8,20,16],
                 [5,64,32],
                 [16,64,5],
                 [5,64,32],
                 [16,64,5],
                 [16,5,32]]
        n_vector = [Float64[],
                    []]
        t_vector = [Float64[],
                    []]
        v = [0.0,0.0,0.0]
        T = 500.0
        H = MaterialLib.enth(MaterialLib.Sodium(), 500.0, p0)

        for ii in 1:NUM_POOL
            pool = SPHParticles.ParticlePool(ii, TAG[ii], NUM[ii])
            cordins = GeoLib.regular(CORDIN[ii][1], CORDIN[ii][2], CORDIN[ii][3], NDIMS[ii])
            SPHParticles.set_var(pool, :rho, rho0)
            SPHParticles.set_var(pool, :V, 1.0/volume)
            SPHParticles.set_var(pool, :m, rho0*volume)
            SPHParticles.set_var(pool, :h, hdx*dx)
            SPHParticles.set_var(pool, :p, p0)
            SPHParticles.set_var(pool, :T, T[ii])
            SPHParticles.set_var(pool, :H, H[ii])
            SPHParticles.set_var(pool, :vx, v[1])
            SPHParticles.set_var(pool, :vy, v[2])
            SPHParticles.set_var(pool, :vz, v[3])
            SPHParticles.set_var(pool, :x, cordins[1])
            SPHParticles.set_var(pool, :y, cordins[2])
            SPHParticles.set_var(pool, :z, cordins[3])
            if TAG[ii] == :Solid
                r = sqrt.(cordins[1].*cordins[1] + cordins[2].*cordins[2])
                SPHParticles.set_var(pool, :nx, -cordins[1]./r)
                SPHParticles.set_var(pool, :ny, -cordins[2]./r)
                SPHParticles.set_var(pool, :nz, 0.0)
                SPHParticles.set_var(pool, :tx, -cordins[2]./r)
                SPHParticles.set_var(pool, :ty, cordins[1]./r)
                SPHParticles.set_var(pool, :tz, 0.0)
            end
            if TAG[ii] == :Inlet
                SPHParticles.set_var(pool, :x, cordins[1].+(dt*v[ii][1]))
                SPHParticles.set_var(pool, :y, cordins[2].+(dt*v[ii][2]))
                SPHParticles.set_var(pool, :z, cordins[3].+(dt*v[ii][3]))
            end
            SPHParticles.add(p, TAG[ii], pool)
        end

        # Set Boundarys
        bounds_1 = GeoLib.Boundary[]
        push!(bounds_1, GeoLib.Boundary(GeoLib.CylinderZ([0.0,0.0], 0.05), <))
        SPHParticles.link([[p[:Solid,1]], [p[:Fluid,1]]], bounds_1)

        bounds_2 = GeoLib.Boundary[]
        push!(bounds_2, GeoLib.Boundary(GeoLib.PlateZ(0.0), <))
        push!(bounds_2, GeoLib.Boundary(GeoLib.PlateZ(0.2), <))
        SPHParticles.link([[p[:Inlet,1]], [p[:Fluid,1]], [p[:Outlet,1]]], bounds_2)
        return p
    end

case = SPHApplication.Case("SPH-TANK", "./cases/output/SPH-tank")
SPHApplication.initialize(case)
SPHApplication.apprun(case)
