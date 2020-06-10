################################################################################
## Heat Conduction Equations
################################################################################
    mutable struct HeatConductionEquationT{DstMat <: MaterialLib.Material, SrcMat <: MaterialLib.Material} <: Equation
        dst::Symbol
        src::Symbol
        dst_mat::DstMat
        src_mat::SrcMat
        delta::Float64
        HeatConductionEquationT(dst::Symbol, src::Symbol, dst_mat=MaterialLib.Sodium(), src_mat=MaterialLib.Sodium(), delta::Float64=0.0) =
        new{typeof(dst_mat), typeof(src_mat)}(dst,src,dst_mat,src_mat,delta)

        HeatConductionEquationT(dstsrc::Symbol, dstsrc_mat, delta::Float64=0.0) =
        HeatConductionEquationT(dstsrc,dstsrc,dstsrc_mat,dstsrc_mat,delta)
    end

    function initialize(eq::HeatConductionEquationT, p::SPHParticles.ParticlePool)
        for d_idx in p.pidx
            p[:aT, d_idx] = 0.0
        end
        #p[:aT] = 0.0
    end

    function loop(eq::HeatConductionEquationT, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool,
                  dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(2)
        DWIJ::Vector{Float64} = zeros(2)
        tmp_mat = zeros(2,2)

        lamdai::Float64 = MaterialLib.lamda(eq.dst_mat, dst[:T, d_idx], dst[:p, d_idx])
        cpi::Float64 = MaterialLib.c_p(eq.dst_mat, dst[:T, d_idx], dst[:p, d_idx])

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
            dst[:aT, d_idx] += kij*delta_T*dot(DWIJ,XIJ)/RIJ/RIJ/(1.0+eq.delta^2)/dst[:rho, d_idx]/src[:V, s_idx]/cpi
        end
    end

    """function loop(eq::HeatConductionEquationT3D, k::SPHKernels.Kernel, src::SPHParticles.ParticlePool,
                  dst::SPHParticles.ParticlePool, d_idx::Int)
        XIJ::Vector{Float64} = zeros(2)
        DWIJ::Vector{Float64} = zeros(2)
        VIJ::Vector{Float64} = zeros(2)
        tmp_mat = zeros(2,2)

        lamdai::Float64 = MaterialLib.lamda(eq.dst_mat, dst[:T, d_idx], dst[:p, d_idx])
        cpi::Float64 = MaterialLib.c_p(eq.dst_mat, dst[:T, d_idx], dst[:p, d_idx])

        for s_idx in src.pidx
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            RIJ = norm(XIJ)
            if RIJ <= 1e-10 continue end
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            tmp_mat += kron(XIJ[1:eq.dim], XIJ[1:2]')*norm(DWIJ)/RIJ/src[:V, s_idx]
        end
        Ra = inv(tmp_mat)

        for s_idx in src.pidx
            HIJ::Float64 = 0.5*(dst[:h, d_idx] + src[:h, s_idx])
            XIJ[1] = dst[:x, d_idx] - src[:x, s_idx]
            XIJ[2] = dst[:y, d_idx] - src[:y, s_idx]
            XIJ[3] = dst[:z, d_idx] - src[:z, s_idx]
            RIJ = norm(XIJ)
            if RIJ <= 1e-10 continue end
            SPHKernels.gradient(k, RIJ, HIJ, XIJ, DWIJ)
            # correction
            DWIJ[1:eq.dim] = Ra*DWIJ[1:eq.dim]
            #
            rijdotdwij = dot(DWIJ,XIJ)
            delta_T = dst[:T, d_idx] - src[:T, s_idx]
            lamdaj::Float64 = MaterialLib.lamda(eq.src_mat, src[:T, s_idx], src[:p, s_idx])
            kij = 4.0*lamdai*lamdaj/(lamdai+lamdaj)
            dst[:aT, d_idx] += src[:m, s_idx]*kij*delta_T*rijdotdwij/R2IJ/(1.0+eq.delta^2)/dst[:rho, d_idx]/src[:rho, s_idx]/cpi
            #println(d_idx,"  ",(tmp/cpi),"  ",src[:rho, s_idx],"  ",lamdai," ",kij,"  ",delta_T," ",rijdotdwij," ",(R2IJ+eq.delta^2)," ",dst[:rho, d_idx])
        end
        #println(src.tag,"  ",d_idx,"  ",ra_correction,"  ",dst[:aT, d_idx],"  ",length(src.pidx))
    end"""
