################################################################################
## Energy Equations
################################################################################
    mutable struct EnergyEquationH{DstMat <: MaterialLib.Material, SrcMat <: MaterialLib.Material} <: Equation
        dst::Symbol
        src::Symbol
        dst_mat::DstMat
        src_mat::SrcMat
        delta::Float64
        EnergyEquationH(dst::Symbol, src::Symbol, dst_mat=MaterialLib.Sodium(), src_mat=MaterialLib.Sodium(), delta::Float64=0.0) =
        new{typeof(dst_mat), typeof(src_mat)}(dst,src,dst_mat,src_mat,delta)

        EnergyEquationH(dstsrc::Symbol, dstsrc_mat, delta::Float64=0.0) =
        EnergyEquationH(dstsrc,dstsrc,dstsrc_mat,dstsrc_mat,delta)
    end

    function initialize(eq::EnergyEquationH, p::SPHParticles.ParticlePool)
        for d_idx in p.pidx
            p[:aH, d_idx] = 0.0
        end
        #p[:aT] = 0.0
    end

    function loop(eq::EnergyEquationH, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool,
                  dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(2)
        DWIJ::Vector{Float64} = zeros(2)
        tmp_mat = zeros(2,2)

        dst[:T, d_idx] = MaterialLib.temperature(eq.dst_mat, dst[:H, d_idx], dst[:p, d_idx])
        lamdai::Float64 = MaterialLib.lamda(eq.dst_mat, dst[:T, d_idx], dst[:p, d_idx])

        for s_idx in src.pidx
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            RIJ = norm(XIJ)
            if RIJ <= 1e-10 continue end
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            tmp_mat += kron(XIJ, XIJ')*norm(DWIJ)/RIJ/src[:V, s_idx]
        end
        Ra = inv(tmp_mat)

        for s_idx in src.pidx
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            RIJ = norm(XIJ)
            if RIJ <= 1e-10 continue end
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            # correction
            DWIJ = Ra*DWIJ
            #
            delta_T = dst[:T, d_idx] - src[:T, s_idx]
            lamdaj::Float64 = MaterialLib.lamda(eq.src_mat, src[:T, s_idx], src[:p, s_idx])
            kij = 4.0*lamdai*lamdaj/(lamdai+lamdaj)
            dst[:aH, d_idx] += kij*delta_T*dot(DWIJ,XIJ)/RIJ/RIJ/(1.0+eq.delta^2)/dst[:rho, d_idx]/src[:V, s_idx]
        end
    end
