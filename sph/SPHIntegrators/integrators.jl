################################################################################
## Integrators
################################################################################
    abstract type Integrator end
    function timestep(T::Integrator, dt::Float64) end

###############################################################################
# `EulerIntegrator` struct
###############################################################################
    mutable struct EulerIntegrator <: Integrator
        target::Symbol
        stepper::IntegratorStep
        particles::SPHParticles.ParticlePoolUnion
        eqs::SPHEquations.Group
        kernel::SPHKernels.Kernel
        EulerIntegrator(tags=:Fluid, stepper=EulerStep()) = new(tags, stepper)
    end

    function timestep(int::EulerIntegrator, dt::Float64)
        SPHEquations.calculate_force(int.eqs, int.particles, int.kernel)
        @time foreach(pool->stage1(int.stepper, pool, dt), int.particles[int.target])
    end

###############################################################################
# `PECIntegrator` struct
###############################################################################
    "In the Predict-Evaluate-Correct (PEC) mode, the system is advanced using

    .. math

    y^{n+frac{1}{2}} = y^n + frac{Delta t}{2}F(y^{n-frac{1}{2}})
    --> Predict

    F(y^{n+frac{1}{2}}) --> Evaluate

    y^{n + 1} = y^n + Delta t F(y^{n+frac{1}{2}})
    "

    mutable struct PECIntegrator <: Integrator
        target::Symbol
        stepper::IntegratorStep
        particles::SPHParticles.ParticlePoolUnion
        eqs::SPHEquations.Group
        kernel::SPHKernels.Kernel
        PECIntegrator(tag=:Fluid, stepper=WCSPHStep()) = new(tag, stepper)
    end

    function timestep(int::PECIntegrator, dt::Float64)
        foreach(pool->initialize(int.stepper, pool), int.particles[int.target])
        # Predict
        foreach(pool->stage1(int.stepper, pool, dt), int.particles[int.target])
        # Call any post-stage functions.
        #do_post_stage(pool, 0.5*dt, 1)
        # Correct
        SPHEquations.calculate_force(int.eqs, int.particles, int.kernel)
        foreach(pool->stage2(int.stepper, pool, dt), int.particles[int.target])
        # Call any post-stage functions.
        #do_post_stage(p, dt, 2)
    end
