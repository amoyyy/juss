## Equations of Boundary
################################################################################
    struct SetWallVelocity <: Equation
        dst::Symbol
        src::Symbol
        SetWallVelocity(dst::Symbol, src::Symbol) = new(dst, src)
        SetWallVelocity() = SetWallVelocity(:Solid, :Fluid)
    end

    function initialize(eq::SetWallVelocity, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:vxf, idx] = 0.0
            p[:vyf, idx] = 0.0
            p[:vzf, idx] = 0.0
            p[:wij, idx] = 0.0
        end
    end

    function loop(eq::SetWallVelocity, k::SPHKernels.Kernel,
                  src::SPHParticles.ParticlePool, dst::SPHParticles.ParticlePool, d_idx::Int)
        # normalisation factor is different from 'V' as the particles
        # near the boundary do not have full kernel support
        XIJ::Vector{Float64} = zeros(3)
        for s_idx in src.pidx
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = sqrt(R2IJ)
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            wij = SPHKernels.kernel(k, RIJ, HIJ)
            dst[:wij, d_idx] += wij
            # sum in Eq. (22)
            # this will be normalized in post loop
            dst[:vxf, d_idx] += src[:vx, s_idx] * wij
            dst[:vyf, d_idx] += src[:vy, s_idx] * wij
            dst[:vzf, d_idx] += src[:vz, s_idx] * wij
        end
    end

    function post_loop(eq::SetWallVelocity, p::SPHParticles.ParticlePool)
        # calculation is done only for the relevant boundary particles.
        # d_wij (and d_uf) is 0 for particles sufficiently away from the
        # solid-fluid interface
        for d_idx in p.pidx
            if p[:wij, d_idx] > 1e-12
                p[:vxf, d_idx] /= p[:wij, d_idx]
                p[:vyf, d_idx] /= p[:wij, d_idx]
                p[:vzf, d_idx] /= p[:wij, d_idx]
            end

            # Dummy velocities at the ghost points using Eq. (23),
            # d_u, d_v, d_w are the prescribed wall velocities.
            p[:vxg, d_idx] = 2*p[:vx, d_idx] - p[:vxf, d_idx]
            p[:vyg, d_idx] = 2*p[:vy, d_idx] - p[:vyf, d_idx]
            p[:vzg, d_idx] = 2*p[:vz, d_idx] - p[:vzf, d_idx]
        end
    end

    #############################################################################
    struct WallNoSlipBC <: Equation
        """
        Parameters
        ----------
        nu : float
            kinematic viscosity

        Notes
        -----
        For this equation the destination particle array should be the
        fluid and the source should be ghost or boundary particles. The
        boundary particles must functionine a prescribed velocity :math:`u_0,
        v_0, w_0`
        """
        dst::Symbol
        src::Symbol
        nu::Float64
        WallNoSlipBC(dst::Symbol, src::Symbol, nu::Float64) = new(dst,src,nu)
        WallNoSlipBC(nu::Float64) = WallNoSlipBC(:Fluid, :Solid, nu)
    end

    function initialize(eq::WallNoSlipBC, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:ax, idx] = 0.0
            p[:ay, idx] = 0.0
            p[:az, idx] = 0.0
        end
    end

    function loop(eq::WallNoSlipBC, k::SPHKernels.Kernel,
                  src::SPHParticles.ParticlePool, dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        DWIJ::Vector{Float64} = zeros(3)
"""        dwijkronrij_sum = zeros(2,2)
        ra_correction = zeros(2,2)
        for s_idx in src.pidx
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = sqrt(R2IJ)
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            dwijkronrij_sum += kron(DWIJ[1:2], XIJ[1:2]')*src[:m, s_idx]/src[:rho, s_idx]
        end
        ra_correction = inv(dwijkronrij_sum)"""
        for s_idx in src.pidx
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = sqrt(R2IJ)
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            EPS = 0.0#1*HIJ*HIJ
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
            #DWIJ[1:2] = ra_correction*DWIJ[1:2]
            Fij = XIJ[1]*DWIJ[1] + XIJ[2]*DWIJ[2] + XIJ[3]*DWIJ[3]

            # viscous contribution (third term) from Eq. (8), with VIJ
            # functionined appropriately using the ghost values
            tmp = 1.0/dst[:m, d_idx] * (Vi2 + Vj2) * (etaij * Fij/(R2IJ + EPS))

            dst[:ax, d_idx] += tmp * (dst[:vx, d_idx] - src[:vxg, s_idx])
            dst[:ay, d_idx] += tmp * (dst[:vy, d_idx] - src[:vyg, s_idx])
            dst[:az, d_idx] += tmp * (dst[:vz, d_idx] - src[:vzg, s_idx])
            #println(d_idx,"  ",tmp,"  ",dst[:vx, d_idx]," ",src[:vxg, s_idx],"  ",dst[:ax, d_idx])
            #println(d_idx,"  ",src[:vxg, s_idx],"  ",dst[:ax, d_idx])
        end
    end

    #############################################################################
    mutable struct SolidWallPressureBC <: Equation
"""        Parameters
        ----------
        rho0 : float
            reference density
        p0 : float
            reference pressure
        b : float
            constant (functionault 1.0)
        gx : float
            Body force per unit mass along the x-axis
        gy : float
            Body force per unit mass along the y-axis
        gz : float
            Body force per unit mass along the z-axis"""
        dst::Symbol
        src::Symbol
        rho0::Float64
        p0::Float64
        b::Float64
        gx::Float64
        gy::Float64
        gz::Float64
        SolidWallPressureBC(dst::Symbol, src::Symbol, rho0::Float64, p0::Float64, gx::Float64, gy::Float64, gz::Float64) =
        new(dst, src, rho0, p0, 1.0, gx, gy, gz)
        SolidWallPressureBC(rho0::Float64, p0::Float64, gx::Float64, gy::Float64, gz::Float64) =
        SolidWallPressureBC(:Solid, :Fluid, rho0, p0, gx, gy, gz)
    end

    function initialize(eq::SolidWallPressureBC, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:p, idx] = 0.0
            p[:wij, idx] = 0.0
        end
    end

    function loop(eq::SolidWallPressureBC, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool, dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        for s_idx in src.pidx
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = sqrt(R2IJ)
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            wij = SPHKernels.kernel(k, RIJ, HIJ)*src[:m, s_idx]/src[:rho, s_idx]
            # numerator of Eq. (27) ax, ay and az are the prescribed wall
            # accelerations which must be functionined for the wall boundary
            # particle
            gdotxij = (eq.gx-dst[:ax, d_idx])*XIJ[1]+(eq.gy-dst[:ay, d_idx])*XIJ[2]+(eq.gz-dst[:az, d_idx])*XIJ[3]
            dst[:p, d_idx] += (src[:p, s_idx]*wij + src[:rho, s_idx]*gdotxij*wij)
            # denominator of Eq. (27)
            dst[:wij, d_idx] += wij
            #if gdotxij != 0.0
            #    println("solid: ",d_idx,"  ",s_idx,"  ",dst[:p, d_idx],"  ",src[:p, s_idx],"  ",wij,"  ",gdotxij)
            #end
        end
    end

    function post_loop(eq::SolidWallPressureBC, p::SPHParticles.ParticlePool)
        for d_idx in p.pidx
            # extrapolated pressure at the ghost particle
            if p[:wij, d_idx] > 1e-14
                p[:p, d_idx] /= p[:wij, d_idx]
            end
            # update the density from the pressure Eq. (28)
            p[:rho, d_idx] = eq.rho0 * (p[:p, d_idx]/eq.p0 + eq.b)
            #if p[:p, d_idx] != 0.0
            #    println("solid: ",d_idx,"  ",p[:p, d_idx],"  ",p[:rho, d_idx])
            #end
        end
    end

    #############################################################################
    struct MonaghanBoundaryForce <: Equation
        dst::Symbol
        src::Symbol
        d_ap::Float64
        MonaghanBoundaryForce(dst::Symbol, src::Symbol, d_ap::Float64=0) =
        new(dst, src, d_ap)
        MonaghanBoundaryForce(d_ap::Float64) = MonaghanBoundaryForce(:Fluid, :Solid, d_ap)
    end

    function loop(eq::MonaghanBoundaryForce, k::SPHKernels.Kernel,
                  src::SPHParticles.ParticlePool, dst::SPHParticles.ParticlePool, d_idx::Int)
        norm::Vector{Float64} = zeros(3)
        tang::Vector{Float64} = zeros(3)
        XIJ::Vector{Float64} = zeros(3)
        for s_idx in src.pidx
            mi = dst[:m, d_idx]
            mj = src[:m, s_idx]
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            # particle sound speed
            cs = dst[:cs, d_idx]

            # boundary normals
            norm[1] = src[:nx, s_idx]
            norm[2] = src[:ny, s_idx]
            norm[3] = src[:nz, s_idx]

            # boundary tangents
            tang[1] = src[:tx, s_idx]
            tang[2] = src[:ty, s_idx]
            tang[3] = src[:tz, s_idx]

            # x and y projections
            x = XIJ[1]*tang[1] + XIJ[2]*tang[2] + XIJ[3]*tang[3]
            y = XIJ[1]*norm[1] + XIJ[2]*norm[2] + XIJ[3]*norm[3]

            # compute the force
            force = 0.0
            q = y/dst[:h, d_idx]
            xabs = abs(x)

            if 0 <= xabs <= eq.d_ap
                beta = 0.02 * cs * cs/y
                tforce = 1.0 - xabs/eq.d_ap
                if 0 < q <= 2.0/3.0
                    nforce =  2.0/3.0
                elseif 2.0/3.0 < q <= 1.0
                    nforce = 2*q*(1.0 - 0.75*q)
                elseif 1.0 < q <= 2.0
                    nforce = 0.5 * (2-q)*(2-q)
                else
                    nforce = 0.0
                end
                force = (mi/(mi+mj)) * nforce * tforce * beta
            end

            # boundary force accelerations
            dst[:ax, d_idx] += force * XIJ[1]
            dst[:ay, d_idx] += force * XIJ[2]
            dst[:az, d_idx] += force * XIJ[3]
        end
    end
    #---------------------------------------------------------------------------
    function wendland_quintic(rij::Float64=1.0, h::Float64=1.0)
        q = rij/h
        q1 = 2.0 - q
        val = 0.0
        # DEBUG
        if 0.0< q < 2.0
            val = (1.0 + 2.5*q + 2*q*q)*q1*q1*q1*q1*q1
        end
        return val
    end

    struct MonaghanKajtarBoundaryForce <: Equation
        dst::Symbol
        src::Symbol
        K::Float64
        beta::Float64
        h::Float64
        MonaghanKajtarBoundaryForce(dst::Symbol=:Fluid, src::Symbol=:Solid, K::Float64=NaN, beta::Float64=NaN, h::Float64=NaN) =
        new(dst, src, K, beta, h)
    end

    function loop(eq::MonaghanKajtarBoundaryForce,
                  d_idx, s_idx,
                  m, ax, ay, az,
                  RIJ, R2IJ, XIJ)
        ma = m[d_idx]
        mb = m[s_idx]
        w = wendland_quintic(RIJ, eq.h)
        force = eq.K/eq.beta * w/R2IJ * 2*mb/(ma + mb)
        ax[d_idx] += force * XIJ[1]
        ay[d_idx] += force * XIJ[2]
        az[d_idx] += force * XIJ[3]
    end

    function loop(eq::MonaghanKajtarBoundaryForce, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool,
                      dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        for s_idx in src.pidx
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = sqrt(R2IJ)
            ma = dst[:m, d_idx]
            mb = src[:m, s_idx]
            w = wendland_quintic(RIJ, eq.h)
            force = eq.K/eq.beta * w/R2IJ * 2*mb/(ma + mb)
            #println(force, XIJ)
            dst[:ax, d_idx] += force * XIJ[1]
            dst[:ay, d_idx] += force * XIJ[2]
            dst[:az, d_idx] += force * XIJ[3]
        end
        #println(d_idx,"  ", dst[:ax, d_idx], "  ",dst[:ay, d_idx],"  ",length(src.pidx))
    end
