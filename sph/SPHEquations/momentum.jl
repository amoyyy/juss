#########################################################################
## MomentumEquationPressureGradient Equations
#########################################################################
    mutable struct MomentumEquation <: Equation
        dst::Symbol
        src::Symbol
        c0::Float64
        alpha::Float64
        beta::Float64
        gx::Float64
        gy::Float64
        gz::Float64
        tensile_correction::Bool
        #MomentumEquation(dst, src, c0=1.0, alpha=1.0, beta=2.0, gx=0.0, gy=0.0, gz=0.0, tc=false) =
        MomentumEquation(dst::Symbol, src::Symbol, c0::Float64, alpha::Float64, beta::Float64, gx::Float64, gy::Float64, gz::Float64, tc::Bool) =
        new(dst, src, c0, alpha, beta, gx, gy, gz, tc)

        MomentumEquation(dstsrc::Symbol, c0::Float64, alpha::Float64, beta::Float64, gx::Float64, gy::Float64, gz::Float64) =
        MomentumEquation(dstsrc, dstsrc, c0, alpha, beta, gx, gy, gz, true)
    end

    """
    MomentumEquation
    Parameters
    --------------------------------------------------------------------------------------------------------------------
    c0  Float64
        reference speed of sound
    alpha : Float64
        produces a shear and bulk viscosity
    beta : Float64
        used to handle high Mach number shocks
    gx : Float64
        body force per unit mass along the x-axis
    gy : Float64
        body force per unit mass along the y-axis
    gz : Float64
        body force per unit mass along the z-axis
    tensilte_correction : bool
        switch for tensile instability correction (functionault: False)
    --------------------------------------------------------------------------------------------------------------------
    """

    function initialize(eq::MomentumEquation, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:ax, idx] = 0.0
            p[:ay, idx] = 0.0
            p[:az, idx] = 0.0
        end
        #p[:dt_cfl] .= 0.0
    end

    # DEBUG
    function loop(eq::MomentumEquation, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool,
                  dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        DWIJ::Vector{Float64} = zeros(3)
        VIJ::Vector{Float64} = zeros(3)

        tmp_mat = zeros(2,2)
        for s_idx in src.pidx
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            RIJ = norm(XIJ)
            if RIJ <= 1e-10 continue end
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            tmp_mat += kron(XIJ[1:2], XIJ[1:2]')*norm(DWIJ)/RIJ/src[:V, s_idx]
        end
        #Ra = inv(tmp_mat)
        for s_idx in src.pidx
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = norm(XIJ)
            RHOIJ = 0.5*(dst[:rho, d_idx] + src[:rho, s_idx])
            RHOIJ1 = 1.0/RHOIJ
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            #
            #print(tmp_mat,"  ",DWIJ,"   ")
            #DWIJ[1:2] = LAPACK.sysv!('U', tmp_mat, DWIJ[1:2])[1]
            #println(DWIJ)
            #DWIJ = Ra*DWIJ
            #
            VIJ[1] = dst[:vx, d_idx] - src[:vx, s_idx]
            VIJ[2] = dst[:vy, d_idx] - src[:vy, s_idx]
            VIJ[3] = dst[:vz, d_idx] - src[:vz, s_idx]
            EPS::Float64 = 0.01 * HIJ * HIJ
            rhoi21 = 1.0/(dst[:rho, d_idx]*dst[:rho, d_idx])
            rhoj21 = 1.0/(src[:rho, s_idx]*src[:rho, s_idx])
            vijdotxij = VIJ[1]*XIJ[1] + VIJ[2]*XIJ[2] + VIJ[3]*XIJ[3]
            piij = 0
            if vijdotxij < 0
                cij = 0.5 * (dst[:cs, d_idx] + src[:cs, s_idx])
                muij = (HIJ * vijdotxij)/(R2IJ + EPS)
                piij = -eq.alpha*cij*muij + eq.beta*muij*muij
                piij = piij*RHOIJ1
            end
            tmpi = dst[:p, d_idx]*rhoi21
            tmpj = src[:p, s_idx]*rhoj21

            ## WDP
            deltap = SPHKernels.get_coef(k)
            WIJ = SPHKernels.kernel(k, RIJ, HIJ)
            WDP = SPHKernels.kernel(k, HIJ*deltap, HIJ)
            fij = WIJ/WDP
            Ri = 0.0
            Rj = 0.0

            # tensile instability correction
            if eq.tensile_correction
                fij = fij*fij
                fij = fij*fij
                if dst[:p, d_idx] > 0
                    Ri = 0.01 * tmpi
                else
                    Ri = 0.2*abs(tmpi)
                end
                if src[:p, s_idx] > 0
                    Rj = 0.01 * tmpj
                else
                    Rj = 0.2 * abs(tmpj)
                end
            end

            # compute the CFL time step factor
            #if R2IJ > 1e-12
            #    _dt_cfl = abs(HIJ * vijdotxij/R2IJ) + eq.c0
            #    dst[:dt_cfl, d_idx] = max(_dt_cfl, dst[:dt_cfl, d_idx])
            #end

            # SPHKernels.gradient and correction terms
            tmp = (tmpi + tmpj) + (Ri + Rj)*fij
            # addup to velocity acceleration
            dst[:ax, d_idx] += -src[:m, s_idx] * (tmp + piij) * DWIJ[1]
            dst[:ay, d_idx] += -src[:m, s_idx] * (tmp + piij) * DWIJ[2]
            dst[:az, d_idx] += -src[:m, s_idx] * (tmp + piij) * DWIJ[3]
        end
    end

    function post_loop(eq::MomentumEquation, p::SPHParticles.ParticlePool)
        p[:ax] .+= eq.gx
        p[:ay] .+= eq.gy
        p[:az] .+= eq.gz

        #acc = sqrt.(p[:ax].*p[:ax] + p[:ay].*p[:ay] + p[:az].*p[:az])

        #p[:dt_force] .= acc
    end

    ##########################################################################################
    mutable struct MomentumEquationPressureGradient <: Equation
        dst::Symbol
        src::Symbol
        pb::Float64
        gx::Float64
        gy::Float64
        gz::Float64
        tdamp::Float64
        MomentumEquationPressureGradient(dst::Symbol, src::Symbol, pb::Float64, gx::Float64, gy::Float64, gz::Float64, tdamp::Float64=0.0) =
        new(dst, src, pb, gx, gy, gz, tdamp)
    end

    function initialize(eq::MomentumEquationPressureGradient, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:ax, idx] = 0.0
            p[:ay, idx] = 0.0
            p[:az, idx] = 0.0

            p[:ax0, idx] = 0.0
            p[:ay0, idx] = 0.0
            p[:az0, idx] = 0.0
        end
    end

    function loop(eq::MomentumEquationPressureGradient, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool, dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        DWIJ::Vector{Float64} = zeros(3)
        for s_idx in src.pidx
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = norm(XIJ)
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            # averaged pressure Eq. (7)
            rhoi = dst[:rho, d_idx]
            rhoj = src[:rho, s_idx]
            p_i = dst[:p, d_idx]
            p_j = src[:p, s_idx]

            pij = rhoj * p_i + rhoi * p_j
            pij /= (rhoj + rhoi)

            # particle volumes; d_V is inverse volume
            Vi = 1.0/dst[:V, d_idx]
            Vj = 1.0/src[:V, s_idx]
            Vi2 = Vi * Vi
            Vj2 = Vj * Vj

            # inverse mass of destination particle
            mi1 = 1.0/dst[:m, d_idx]

            # accelerations 1st term in Eq. (8)
            tmp = -pij * mi1 * (Vi2 + Vj2)

            dst[:ax, d_idx] += tmp * DWIJ[1]
            dst[:ay, d_idx] += tmp * DWIJ[2]
            dst[:az, d_idx] += tmp * DWIJ[3]

            # contribution due to the background pressure Eq. (13)
            tmp = -eq.pb * mi1 * (Vi2 + Vj2)

            dst[:ax0, d_idx] += tmp * DWIJ[1]
            dst[:ay0, d_idx] += tmp * DWIJ[2]
            dst[:az0, d_idx] += tmp * DWIJ[3]
        end
    end

    function post_loop(eq::MomentumEquationPressureGradient, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool, dst::SPHParticles.ParticlePool, d_idx::Int)
        # damped accelerations due to body or external force
        damping_factor = 1.0
        if t < eq.tdamp
            damping_factor = 0.5 * (sin((-0.5 + t/eq.tdamp)*pi) + 1.0)
        end

        dst[:ax, d_idx] += eq.gx * damping_factor
        dst[:ay, d_idx] += eq.gy * damping_factor
        dst[:az, d_idx] += eq.gz * damping_factor
    end

#########################################################################
## ArtificialViscosity Equations
#########################################################################
    mutable struct ArtificialViscosity <: Equation
        dst::Symbol
        src::Symbol
        alpha::Float64
        cs::Float64
        ArtificialViscosity(dst::Symbol, src::Symbol, alpha::Float64, cs::Float64) =
        new(dst, src, alpha, cs)

        ArtificialViscosity(dstsrc::Symbol, alpha::Float64, cs::Float64) =
        ArtificialViscosity(dstsrc, dstsrc, alpha, cs)
    end

    function initialize(eq::ArtificialViscosity, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:ax, idx] = 0.0
            p[:ay, idx] = 0.0
            p[:az, idx] = 0.0
        end
    end

    function loop(eq::ArtificialViscosity, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool, dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        DWIJ::Vector{Float64} = zeros(3)
        VIJ::Vector{Float64} = zeros(3)
        for s_idx in src.pidx
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = norm(XIJ)
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            # averaged pressure Eq. (7)
            VIJ[1] = dst[:vx, d_idx] - src[:vx, s_idx]
            VIJ[2] = dst[:vy, d_idx] - src[:vy, s_idx]
            VIJ[3] = dst[:vz, d_idx] - src[:vz, s_idx]
            vijdotrij = VIJ[1]*XIJ[1] + VIJ[2]*XIJ[2] + VIJ[3]*XIJ[3]
            RHOIJ = 0.5*(dst[:rho, d_idx] + src[:rho, s_idx])
            RHOIJ1 = 1.0/RHOIJ
            # scalar part of the accelerations Eq. (11)
            piij::Float64 = 0.0
            EPS::Float64 = 0.01 * HIJ * HIJ
            if vijdotrij < 0
                muij = (HIJ * vijdotrij)/(R2IJ + EPS)
                piij = -eq.alpha*eq.cs*muij
                piij = src[:m, s_idx]*piij*RHOIJ1
            end
            dst[:ax, d_idx] += -piij * DWIJ[1]
            dst[:ay, d_idx] += -piij * DWIJ[2]
            dst[:az, d_idx] += -piij * DWIJ[3]
        end
    end

#########################################################################
## Viscosity Equations
#########################################################################
    struct Viscosity <: Equation
        dst::Symbol
        src::Symbol
        nu::Float64
        Viscosity(dst::Symbol, src::Symbol, nu::Float64) = new(dst,src,nu)
        Viscosity(nu::Float64) = Viscosity(:Fluid, :Fluid, nu)
    end

    function initialize(eq::Viscosity, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:ax, idx] = 0.0
            p[:ay, idx] = 0.0
            p[:az, idx] = 0.0
        end
    end

    function loop(eq::Viscosity, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool, dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        VIJ::Vector{Float64} = zeros(3)
        DWIJ::Vector{Float64} = zeros(3)
        for s_idx in src.pidx
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = sqrt(R2IJ)
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            EPS = 0.01*HIJ*HIJ
            # averaged shear viscosity Eq. (6).
            etai = eq.nu * dst[:rho, d_idx]
            etaj = eq.nu * src[:rho, s_idx]
            etaij = 2 * (etai * etaj)/(etai + etaj)

            # particle volumes; d_V inverse volume.
            Vi = 1.0/dst[:V, d_idx]
            Vj = 1.0/src[:V, s_idx]
            Vi2 = Vi * Vi
            Vj2 = Vj * Vj

            # scalar part of the kernel gradient
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            Fij = XIJ[1]*DWIJ[1] + XIJ[2]*DWIJ[2] + XIJ[3]*DWIJ[3]

            # viscous contribution (third term) from Eq. (8), with VIJ
            # functionined appropriately using the ghost values
            tmp = 1.0/dst[:m, d_idx] * (Vi2 + Vj2) * (etaij * Fij/(R2IJ + EPS))

            dst[:ax, d_idx] += tmp * (dst[:vx, d_idx] - src[:vx, s_idx])
            dst[:ay, d_idx] += tmp * (dst[:vy, d_idx] - src[:vy, s_idx])
            dst[:az, d_idx] += tmp * (dst[:vz, d_idx] - src[:vz, s_idx])
        end
    end

#########################################################################
## MomentumEquationArtificialStress Equations
#########################################################################
    struct ArtificialStress <: Equation
        dst::Symbol
        src::Symbol
        ArtificialStress(dst::Symbol, src::Symbol) = new(dst, src)
        ArtificialStress() = ArtificialStress(:Fluid, :Fluid)
    end

    function initialize(eq::ArtificialStress, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:ax, idx] = 0.0
            p[:ay, idx] = 0.0
            p[:az, idx] = 0.0
        end
    end

    function loop(eq::ArtificialStress, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool, dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        VIJ::Vector{Float64} = zeros(3)
        DWIJ::Vector{Float64} = zeros(3)
        #Ai::Vector{Float64} = zeros(3,3)
        #Aj::Vector{Float64} = zeros(3,3)

        rhoi = dst[:rho, d_idx]
        # physical and advection velocities
        ui = dst[:vx, d_idx]
        uhati = dst[:vx0, d_idx]
        vi = dst[:vy, d_idx]
        vhati = dst[:vy0, d_idx]
        wi = dst[:vz, d_idx]
        whati = dst[:vz0, d_idx]

        # artificial stress tensor
        Axxi = rhoi*ui*(uhati - ui)
        Axyi = rhoi*ui*(vhati - vi)
        Axzi = rhoi*ui*(whati - wi)
        Ayxi = rhoi*vi*(uhati - ui)
        Ayyi = rhoi*vi*(vhati - vi)
        Ayzi = rhoi*vi*(whati - wi)
        Azxi = rhoi*wi*(uhati - ui)
        Azyi = rhoi*wi*(vhati - vi)
        Azzi = rhoi*wi*(whati - wi)

        for s_idx in src.pidx
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = norm(XIJ)
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)

            rhoj = src[:rho, s_idx]
            uj = src[:vx, s_idx]
            uhatj = src[:vx0, s_idx]
            vj = src[:vy, s_idx]
            vhatj = src[:vy0, s_idx]
            wj = src[:vz, s_idx]
            whatj = src[:vz0, s_idx]

            Axxj = rhoj*uj*(uhatj - uj)
            Axyj = rhoj*uj*(vhatj - vj)
            Axzj = rhoj*uj*(whatj - wj)
            Ayxj = rhoj*vj*(uhatj - uj)
            Ayyj = rhoj*vj*(vhatj - vj)
            Ayzj = rhoj*vj*(whatj - wj)
            Azxj = rhoj*wj*(uhatj - uj)
            Azyj = rhoj*wj*(vhatj - vj)
            Azzj = rhoj*wj*(whatj - wj)

            # particle volumes; d_V is inverse volume.
            Vi = 1.0/dst[:V, d_idx]
            Vj = 1.0/src[:V, s_idx]
            Vi2 = Vi * Vi
            Vj2 = Vj * Vj

            # contraction of stress tensor with kernel gradient
            Ax = 0.5*(
                (Axxi + Axxj)*DWIJ[1] +
                (Axyi + Axyj)*DWIJ[2] +
                (Axzi + Axzj)*DWIJ[3]
            )

            Ay = 0.5*(
                (Ayxi + Ayxj)*DWIJ[1] +
                (Ayyi + Ayyj)*DWIJ[2] +
                (Ayzi + Ayzj)*DWIJ[3]
            )

            Az = 0.5*(
                (Azxi + Azxj)*DWIJ[1] +
                (Azyi + Azyj)*DWIJ[2] +
                (Azzi + Azzj)*DWIJ[3]
            )

            # accelerations 2nd part of Eq. (8)
            tmp = 1.0/dst[:m, d_idx] * (Vi2 + Vj2)

            dst[:ax, d_idx] += tmp * Ax
            dst[:ay, d_idx] += tmp * Ay
            dst[:az, d_idx] += tmp * Az
        end
    end

#########################################################################
## XSPHCorrection Equations
#########################################################################
    struct XSPHCorrection <: Equation
        dst::Symbol
        src::Symbol
        eps::Float64
        XSPHCorrection(dst::Symbol, src::Symbol, eps::Float64) = new(dst, src, eps)
        XSPHCorrection(eps::Float64) = XSPHCorrection(:Fluid, :Fluid, eps)
    end

    function initialize(eq::XSPHCorrection, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:ax, idx] = 0.0
            p[:ay, idx] = 0.0
            p[:az, idx] = 0.0
        end
    end

    function loop(eq::XSPHCorrection, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool, dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        VIJ::Vector{Float64} = zeros(3)
        for s_idx in src.pidx
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = norm(XIJ)
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            WIJ = SPHKernels.kernel(k, RIJ, HIJ)
            RHOIJ = 0.5*(dst[:rho, d_idx] + src[:rho, s_idx])
            RHOIJ1 = 1.0/RHOIJ
            tmp = -eq.eps * src[:m, s_idx]*WIJ*RHOIJ1
            VIJ[1] = dst[:vx, d_idx] - src[:vx, s_idx]
            VIJ[2] = dst[:vy, d_idx] - src[:vy, s_idx]
            VIJ[3] = dst[:vz, d_idx] - src[:vz, s_idx]
            dst[:ax, d_idx] += tmp * VIJ[1]
            dst[:ay, d_idx] += tmp * VIJ[2]
            dst[:az, d_idx] += tmp * VIJ[3]
        end
    end

    function post_loop(eq::XSPHCorrection, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:ax, idx] += p[:vx, idx]
            p[:ay, idx] += p[:vy, idx]
            p[:az, idx] += p[:vz, idx]
        end
    end
