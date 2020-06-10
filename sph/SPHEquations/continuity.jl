################################################################################
## Continuity Equations
################################################################################
    mutable struct NumberDensity <: Equation
        dst::Symbol
        src::Symbol
        NumberDensity(dst, src) = new(dst,src)
        NumberDensity(dstsrc=:Fluid) = NumberDensity(dstsrc,dstsrc)
    end

    function initialize(eq::NumberDensity, p::SPHParticles.ParticlePool) nothing end

    function loop(eq::NumberDensity, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool,
                  dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        wij_correct_coeff::Float64 = 0.0
        rho0::Float64 = 0.0
        V0::Float64  = 0.0
        for s_idx in src.pidx
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = sqrt(R2IJ)
            wij = SPHKernels.kernel(k, RIJ, HIJ)
            wij_correct_coeff += wij*src[:m, s_idx]/src[:rho, s_idx]
            #println("  ",d_idx,"  ",wij,"  ",src[:m, s_idx]," ",src[:rho, s_idx])
            #wij *= src[:m, s_idx]/src[:rho, s_idx]
            rho0 += src[:m, s_idx]*wij
            # V is the inverse volume
            V0 += wij
        end
        dst[:rho, d_idx] = rho0/wij_correct_coeff
        dst[:V, d_idx] = V0/wij_correct_coeff
    end

    ##########################################################################################
    mutable struct VolumeSummation <: Equation
        dst::Symbol
        src::Symbol
        VolumeSummation(dst, src) = new(dst,src)
        VolumeSummation(dstsrc=:Fluid) = NumberDensity(dstsrc,dstsrc)
    end

    function initialize(eq::VolumeSummation, p::SPHParticles.ParticlePool) nothing  end

    function loop(eq::VolumeSummation, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool,
                  dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        wij_correct_coeff::Float64 = 0.0
        V0::Float64  = 0.0
        for s_idx in src.pidx
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = sqrt(R2IJ)
            wij = SPHKernels.kernel(k, RIJ, HIJ)
            wij_correct_coeff += wij*src[:m, s_idx]/src[:rho, s_idx]
            # V is the inverse volume
            V0 += wij
        end
        dst[:V, d_idx] = V0/wij_correct_coeff
    end

    ##########################################################################################
    struct VolumeFromMassDensity <: Equation
        dst::Symbol
        src::Symbol
        VolumeFromMassDensity(dst::Symbol) = new(dst, :NONE)
    end

    function loop(eq::VolumeFromMassDensity, p::SPHParticles.ParticlePool, idx::Int64)
        """**Set the inverse volume using mass density**"""
        #for idx in p.pidx
            p[:V, idx] = p[:rho, idx]/p[:m, idx]
        #end
    end

    ##########################################################################################
    mutable struct ContinuityEquation <: Equation
        dst::Symbol
        src::Symbol
        ContinuityEquation(dst, src) = new(dst,src)
        ContinuityEquation(dstsrc=:Fluid) = ContinuityEquation(dstsrc,dstsrc)
    end

    function initialize(eq::ContinuityEquation, p::SPHParticles.ParticlePool)
        for idx in p.pidx
            p[:arho, idx] = 0.0
        end
    end

    function loop(eq::ContinuityEquation, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool,
                  dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(3)
        DWIJ::Vector{Float64} = zeros(3)
        VIJ::Vector{Float64} = zeros(3)
        #println(d_idx)
        dwijkronrij_sum = zeros(2,2)
        ra_correction = zeros(2,2)
"""        for s_idx in src.pidx
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
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            R2IJ = XIJ[1]*XIJ[1] + XIJ[2]*XIJ[2] + XIJ[3]*XIJ[3]
            if R2IJ == 0.0 continue end
            RIJ = sqrt(R2IJ)
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            #
            #DWIJ[1:2] = ra_correction*DWIJ[1:2]
            #
            VIJ[1] = dst[:vx, d_idx] - src[:vx, s_idx]
            VIJ[2] = dst[:vy, d_idx] - src[:vy, s_idx]
            VIJ[3] = dst[:vz, d_idx] - src[:vz, s_idx]
            vijdotdwij = DWIJ[1]*VIJ[1] + DWIJ[2]*VIJ[2] + DWIJ[3]*VIJ[3]
            dst[:arho, d_idx] += src[:m, s_idx]*vijdotdwij
        end
    end
    ##########################################################################################
