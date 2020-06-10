###############################################################################
## NNP functions
###############################################################################
    function is_particle_near(x0::Float64, y0::Float64, x1::Float64, y1::Float64, distance::Float64)
        sum::Float64 = abs2(x0-x1)
        distance2 = abs2(distance)
        if sum > distance2
            return false
        end
        sum += abs2(y0-y1)
        if sum > distance2
            return false
        end
        if sum == 0.0
            return false
        end
        return true
    end

    function is_particle_near(x0::Float64, y0::Float64, z0::Float64, x1::Float64, y1::Float64, z1::Float64, distance::Float64)
        sum::Float64 = abs2(x0-x1)
        distance2 = abs2(distance)
        if sum > distance2
            return false
        end
        sum += abs2(y0-y1)
        if sum > distance2
            return false
        end
        sum += abs2(z0-z1)
        if sum > distance2
            return false
        end
        if sum == 0.0
            return false
        end
        return true
    end

    function get_adjacency_particles(x::Float64, y::Float64, z::Float64, h::Float64, src::ParticlePool, dst::ParticlePool)
        re::Array{Bool,1} = fill(false, src.num_of_particles)
        map!(idx->is_particle_near(x, y, z, src[:x, idx], src[:y, idx], src[:z, idx], h), re, 1:src.num_of_particles)
        if length(src[re]) > 0
            copy(dst, src, src[re])
        end
    end

    function get_adjacency_particles(s_idx::Int, src::ParticlePool, nnps::ParticlePool)
        # pre clear
        #println("Before ", nnps.num_of_particles, "  ",nnps.num_of_reservation)
        clear(nnps)
        #println("After  ", nnps.num_of_particles, "  ",nnps.num_of_reservation)
        flag1::Bool = (nnps.tag == :ALL)
        # temp vars
        x0::Float64 = src[:x, s_idx]
        y0::Float64 = src[:y, s_idx]
        z0::Float64 = src[:z, s_idx]
        h0::Float64 = src[:h, s_idx]
        for pool in src.links
            #flag2::Bool = (pool.tag != :Inlet && pool.tag != :Outlet)
            if flag1 || pool.tag == nnps.tag
                get_adjacency_particles(x0, y0, z0, h0, pool, nnps)
            end
        end
        if flag1 || src.tag == nnps.tag
            get_adjacency_particles(x0, y0, z0, h0, src, nnps)
        end
    end
