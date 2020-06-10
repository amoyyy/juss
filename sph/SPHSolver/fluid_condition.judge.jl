################################################################################
## FluidConditionJudge
################################################################################
    function cal_timestep(dt::Float64, hdx::Float64, pools::Vector{SPHParticles.ParticlePool})
        dt_new::Float64 = dt
        for pool in pools
            dx2 = (minimum(pool["h"])/hdx)^2/2
            dt_h = minimum(pool["rho"].*pool["cp"]./pool["lamda"])*dx2
            dt_v = minimum(pool["h"]./(pool["cs"].+pool["dt_cfl"]))
            dt_f = minimum(sqrt.(pool["h"]./pool["dt_force"]))
            dt_new = minimum(dt_v, dt_f, dt_h, dt_new)
        end
        return dt_new
    end

    function cal_smoothing_length(dim::Int, hdx::Float64, pools::Vector{SPHParticles.ParticlePool})
        dim1 = 1.0/dim
        """
        if minimum(p.data["V"]) < 0
            println("fking error")
        end
        """
        foreach(pool->(pool["h"] .= pool["V"].^(dim1).*hdx), pools)
    end
