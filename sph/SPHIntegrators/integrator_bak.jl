

###############################################################################
# A leap-frog integrator.
###############################################################################
    struct LeapFrogIntegrator <: Integrator end

    function timestep(int::LeapFrogIntegrator, dt::Float64)
        stage1(p)
        do_post_stage(p, 0.5*dt, 1)
        SPHEquations.calculate_force(int.eqs, int.particles, int.kernel)
        stage2(p)
        do_post_stage(p, dt, 2)
    end

###############################################################################
"A Position-Extended Forest-Ruth-Like integrator [Omeylan2002]_

References
----------
.. [Omeylan2002] I.M. Omelyan, I.M. Mryglod and R. Folk, <Optimized
Forest-Ruth- and Suzuki-like algorithms for integration of motion
in many-body systems>, Computer Physics Communications 146, 188 (2002)
http//arxiv.org/abs/cond-mat/0110585

"
    struct PEFRLIntegrator <: Integrator end

    function timestep(int::PEFRLIntegrator, step::IntegratorStep,  p::SPHParticles.ParticlePool, dt::Float64)
        stage1(p)
        do_post_stage(p, 0.1786178958448091*dt, 1)

        SPHEquations.calculate_force(int.eqs, int.particles, kernel)
        stage2(p)
        do_post_stage(p, 0.1123533131749906*dt, 2)

        SPHEquations.calculate_force(int.eqs, int.particles, kernel)
        stage3(p)
        do_post_stage(p, 0.8876466868250094*dt, 3)

        SPHEquations.calculate_force(int.eqs, int.particles, kernel)
        stage4(p)
        do_post_stage(p, 0.8213821041551909*dt, 4)

        SPHEquations.calculate_force(int.eqs, int.particles, kernel)
        stage5(p)
        do_post_stage(p, dt, 5)
    end

###############################################################################
# `EPECIntegrator` struct
###############################################################################
"
Predictor corrector integrators can have two modes of
operation.

In the Evaluate-Predict-Evaluate-Correct (EPEC) mode, the
system is advanced using

.. math

F(y^n) --> Evaluate

y^{n+frac{1}{2}} = y^n + F(y^n) --> Predict

F(y^{n+frac{1}{2}}) --> Evaluate

y^{n+1} = y^n + Delta t F(y^{n+frac{1}{2}}) --> Correct

Notes

The Evaluate stage of the integrator forces a function
evaluation. Therefore, the PEC mode is much faster but relies on
old accelertions for the Prediction stage.

In the EPEC mode, the final corrector can be modified to

math`y^{n+1} = y^n + frac{Delta t}{2}left( F(y^n) +
F(y^{n+frac{1}{2}}) right)`

This would require additional storage for the accelerations.

"
    struct EPECIntegrator <: Integrator end

    function timestep(int::EPECIntegrator, step::IntegratorStep,  p::SPHParticles.ParticlePool, dt::Float64)
        initialize(p)

        SPHEquations.calculate_force(int.eqs, int.particles, kernel)

        # Predict
        stage1(p)

        # Call any post-stage functions.
        do_post_stage(p, 0.5*dt, 1)

        SPHEquations.calculate_force(int.eqs, int.particles, kernel)

        # Correct
        stage2(p)

        # Call any post-stage functions.
        do_post_stage(p, dt, 2)
    end

###############################################################################
# `TVDRK3Integrator` struct
###############################################################################
    "
    In the TVD-RK3 integrator, the system is advanced using

    .. math

    y^{n + frac{1}{3}} = y^n + Delta t F( y^n )

    y^{n + frac{2}{3}} = frac{3}{4}y^n +
    frac{1}{4}(y^{n + frac{1}{3}} + Delta t F(y^{n + frac{1}{3}}))

    y^{n + 1} = frac{1}{3}y^n + frac{2}{3}(y^{n + frac{2}{3}}
    + Delta t F(y^{n + frac{2}{3}}))

    "
    struct TVDRK3Integrator <: Integrator end

    function timestep(int::TVDRK3Integrator, step::IntegratorStep,  p::SPHParticles.ParticlePool, dt::Float64)
        initialize(p)

        # stage 1
        SPHEquations.calculate_force(int.eqs, int.particles, kernel)
        stage1(p)
        do_post_stage(p, 1.0/3*dt, 1)

        # stage 2
        SPHEquations.calculate_force(int.eqs, int.particles, kernel)
        stage2(p)
        do_post_stage(p, 2.0/3*dt, 2)

        # stage 3 and end
        SPHEquations.calculate_force(int.eqs, int.particles, kernel)
        stage3(p)
        do_post_stage(p, dt, 3)
    end

###############################################################################
