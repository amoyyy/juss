    ##--------------------------------------------------------------------------
    ## Equations of State
    mutable struct TaitEOS <: Equation
        """
        Parameters
        ----------
        rho0 : float
            reference density of fluid particles
        c0 : float
            maximum speed of sound expected in the system
        gamma : float
            constant
        p0 : float
            reference pressure in the system (defaults to zero).

        Notes
        -----
        The reference speed of sound, c0, is to be taken approximately as
        10 times the maximum expected velocity in the system. The particle
        sound speed is given by the usual expression:
        """
        dst::Symbol
        src::Symbol
        rho0::Float64
        c0::Float64
        p0::Float64
        gamma::Float64
        rho01::Float64
        gamma1::Float64
        B::Float64
        TaitEOS(dst::Symbol=:Fluid, src::Symbol=:NONE, rho0::Float64=1.0, c0::Float64=100.0, p0::Float64=0.0, gamma::Float64=7.0) =
        new(dst,:NONE, rho0, c0, p0, gamma, 1.0/rho0, 0.5*(gamma-1.0), rho0*c0*c0/gamma)
    end

"""    function initialize(eq::TaitEOS, dst::SPHParticles.ParticlePool)
        for d_idx in dst.pidx
            ratio = dst[:rho, d_idx] * eq.rho01
            tmp = ratio^(eq.gamma)
            dst[:p, d_idx] = (eq.p0 + eq.B * (tmp - 1.0))
        end
    end"""

    function loop(eq::TaitEOS, dst::SPHParticles.ParticlePool, d_idx::Int64)
        ratio = dst[:rho, d_idx] * eq.rho01
        tmp = ratio^(eq.gamma)
        dst[:p, d_idx] = (eq.p0 + eq.B * (tmp - 1.0))
    end
