###############################################################################
# IntegratorSteps
###############################################################################
    abstract type IntegratorStep end
"""    function stage1(T::IntegratorStep) end
    function stage2(T::IntegratorStep) end
    function stage3(T::IntegratorStep) end
    function stage4(T::IntegratorStep) end
    function stage5(T::IntegratorStep) end"""

###############################################################################
# `EulerStep` struct
###############################################################################
    struct EulerStep <: IntegratorStep end
    """Fast but inaccurate integrator. Use this for testing"""
    #function stage1(p::EulerStep, d_idx, d_u, d_v, d_w, d_au, d_av, d_aw, d_x, d_y,
    #              d_z, d_rho, d_arho, dt)
    #    d_u[d_idx] += dt*d_au[d_idx]
    #    d_v[d_idx] += dt*d_av[d_idx]
    #    d_w[d_idx] += dt*d_aw[d_idx]
    #    d_x[d_idx] += dt*d_u[d_idx]
    #    d_y[d_idx] += dt*d_v[d_idx]
    #    d_z[d_idx] += dt*d_w[d_idx]
    #    d_rho[d_idx] += dt*d_arho[d_idx]
    #end
    function stage1(p::EulerStep, pool::SPHParticles.ParticlePool, dt::Float64)
        pool[:vx] += dt*pool[:ax]
        pool[:vy] += dt*pool[:ay]
        pool[:vz] += dt*pool[:az]

        pool[:x] += dt*pool[:vx]
        pool[:y] += dt*pool[:vy]
        pool[:z] += dt*pool[:vz]

    #    pool[:H] += dt*pool[:aH]
        pool[:H] += dt*pool[:aH]
        #println("$(length(d_idx)),  entered integrator step")
    end

###############################################################################
# `WCSPHStep` struct
###############################################################################
    struct WCSPHStep <: IntegratorStep end
    """Standard Predictor Corrector integrator for the WCSPH formulation

    Use this integrator for WCSPH formulations. In the predictor step,
    the particles are advanced to `t + dt/2`. The particles are then
    advanced with the new force computed at this position.

    This integrator can be used in PEC or EPEC mode.

    The same integrator can be used for other problems. Like for
    example solid mechanics (see SolidMechStep)

    """
    function initialize(p::WCSPHStep, pool::SPHParticles.ParticlePool)
        for idx in pool.pidx
            pool[:x0, idx] = pool[:x, idx]
            pool[:y0, idx] = pool[:y, idx]
            pool[:z0, idx] = pool[:z, idx]

            pool[:vx0, idx] = pool[:vx, idx]
            pool[:vy0, idx] = pool[:vy, idx]
            pool[:vz0, idx] = pool[:vz, idx]

            pool[:rho0, idx] = pool[:rho, idx]
        end
    end

    function stage1(p::WCSPHStep, pool::SPHParticles.ParticlePool, dt::Float64)
        dtb2 = 0.5*dt
        pool[:vx] = pool[:vx0] + dtb2*pool[:ax]
        pool[:vy] = pool[:vy0] + dtb2*pool[:ay]
        #pool[:vz] = pool[:vz0] + dtb2*pool[:az]

        pool[:x] = pool[:x0] + dtb2*pool[:vx]
        pool[:y] = pool[:y0] + dtb2*pool[:vy]
        #pool[:z] = pool[:z0] + dtb2*pool[:vz]

        # Update densities and smoothing lengths from the accelerations
        pool[:rho] = pool[:rho0] + dtb2*pool[:arho]
        #pool[:V] = pool[:rho]./pool[:m]
    end

    function stage2(p::WCSPHStep, pool::SPHParticles.ParticlePool, dt::Float64)
        pool[:vx] = pool[:vx0] + dt*pool[:ax]
        pool[:vy] = pool[:vy0] + dt*pool[:ay]
        #pool[:vz] = pool[:vz0] + dt*pool[:az]

        pool[:x] = pool[:x0] + dt*pool[:vx]
        pool[:y] = pool[:y0] + dt*pool[:vy]
        #pool[:z] = pool[:z0] + dt*pool[:vz]

        # Update densities and smoothing lengths from the accelerations
        pool[:rho] = pool[:rho0] + dt*pool[:arho]
        #pool[:V] = pool[:rho]./pool[:m]
    end

###############################################################################
# `TransportVelocityStep` struct
###############################################################################
    struct TransportVelocityStep <: IntegratorStep end
    """Integrator defined in 'A transport velocity formulation for
    smoothed particle hydrodynamics', 2013, JCP, 241, pp 292--307

    For a predictor-corrector style of integrator, this integrator
    should operate only in PEC mode.

    """
    function initialize(p::TransportVelocityStep, pool::SPHParticles.ParticlePool)
        nothing
    end

    function stage1(p::TransportVelocityStep, pool::SPHParticles.ParticlePool, dt::Float64)
        dtb2 = 0.5*dt
        # velocity update eqn (14)
        pool[:vx] += pool[:ax]*dtb2
        pool[:vy] += pool[:ay]*dtb2
        pool[:vz] += pool[:az]*dtb2
        # advection velocity update eqn (15)
        pool[:vx0] += pool[:vx] + pool[:ax0]*dtb2
        pool[:vy0] += pool[:vy] + pool[:ay0]*dtb2
        pool[:vz0] += pool[:vz] + pool[:az0]*dtb2
        # position update eqn (16)
        pool[:x] += pool[:vx0]*dt
        pool[:y] += pool[:vy0]*dt
        pool[:z] += pool[:vz0]*dt
    end

    function stage2(p::TransportVelocityStep, pool::SPHParticles.ParticlePool, dt::Float64)
        dtb2 = 0.5*dt
        # corrector update eqn (17)
        pool[:vx] += pool[:ax]*dtb2
        pool[:vy] += pool[:ay]*dtb2
        pool[:vz] += pool[:az]*dtb2
        # magnitude of velocity squared
        #d_vmag2[d_idx] = (d_u[d_idx]*d_u[d_idx] + d_v[d_idx]*d_v[d_idx] +
        #                  d_w[d_idx]*d_w[d_idx])
    end
