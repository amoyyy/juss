###############################################################################
## SPHParticles/ParticlePoolUnion
###############################################################################
    mutable struct ParticlePoolUnion
        num_of_pools::Int64
        pools::Dict{Symbol,Vector{ParticlePool}}
        comments::String

        # Constructors
        ParticlePoolUnion() =
        new(0, Dict{Symbol,Vector{ParticlePool}}(), "A BRAND NEW SYSTEM!")
    end

##------------------------------------------------------------------------------
## PoolUnion Function
##------------------------------------------------------------------------------
    function Base.getindex(p::ParticlePoolUnion, tag::Symbol)
        return p.pools[tag]
    end

    function Base.getindex(p::ParticlePoolUnion, tag::Symbol, idx::Int)
        return p.pools[tag][idx]
    end

    function Base.keys(p::ParticlePoolUnion)
        return keys(p.pools)
    end

    # low efficiency
    function Base.getindex(p::ParticlePoolUnion)
        return union(values(p.pools))
    end

    ## Modification Functions
    function set_var(p::ParticlePoolUnion, tag::Symbol, var_tags::Vector{String})
        foreach(pool->set_var(pool, var_tags), p[tag])
    end

    function add(p::ParticlePoolUnion, tag::Symbol, pool::ParticlePool)
        @assert pool.tag == tag "ERROR: POOL TAG DISMATCH!"
        if tag in keys(p)
            push!(p[pool.tag], pool)
        else
            push!(p.pools, (pool.tag => [pool]))
        end
        p.num_of_pools += 1
    end

    function clear(p::ParticlePoolUnion)
        for tag in keys(p)
            foreach(clear, p[tag])
        end
    end

    # update pools with specific tag
    function update(p::ParticlePoolUnion, tag::Symbol)
        #if tag in tags
            foreach(update, p.pools[tag])
        #end
    end

    # update all pools
    function update(p::ParticlePoolUnion)
        foreach(tag->foreach(update, p.pools[tag]), keys(p))
    end

    function check(p::ParticlePoolUnion)
        @assert p.num_of_pools==length(p.tags)==length(keys(p)) "ERROR: PoolDict: $(d.id)!"
        foreach(tag->foreach(check, p.pools[tag]), keys(p))
    end

    function Base.display(p::ParticlePoolUnion)
        println("The Number of Pools : $(p.num_of_pools)")
        foreach(tag->foreach(display, p.pools[tag]), keys(p))
    end
