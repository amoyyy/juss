###############################################################################
## SPHParticles/ParticlePool
###############################################################################
    mutable struct ParticlePool
        id::Int
        tag::Symbol
        num_of_particles::Int
        num_of_reservation::Int
        data::Dict{Symbol,Vector{Float64}}
        pidx::Vector{Int}
        ridx::Vector{Int}       # make ridx transparent
        links::Vector{ParticlePool}
        boundarys::Vector{GeoLib.Boundary}

        # Constructors
        ParticlePool(id::Int, tag::Symbol, nop::Int, nor::Int) =
        init(new(id, tag, nop, nor, Dict{Symbol,Vector{Float64}}(),
                 Array(1:nop), Array(nop+1:nop+nor),
                 Vector{ParticlePool}(), Vector{GeoLib.Boundary}()))

        ParticlePool(id::Int, tag::Symbol, nop::Int) = ParticlePool(id, tag, nop, Int(ceil(nop*0.2)))

        # This is the constructor for temporary pool object
        ParticlePool(tag::Symbol) = ParticlePool(-1, tag, 0, 50)
    end

    function init(p::ParticlePool)
        set_var(p, Main.VARSALL)
        #if p.tag == :Solid
        #    set_var(p, union(Main.VARS,Main.VARS_BOUNDARY_ADDITION))
        #elseif p.tag == :Fluid
        #        set_var(p, union(Main.VARS,Main.VARS_FLUID_ADDITION))
        #else
        #    set_var(p, union(Main.VARS,Main.VARS_FLUID_ADDITION))
#            set_var(p, Main.VARS)
        #end
        return p
    end

    ##------------------------------------------------------------------------------
    ## Pools Function
    ##------------------------------------------------------------------------------
    ## Override System Functions
    function Base.getindex(p::ParticlePool, tag::Symbol)
        return p.data[tag][p.pidx]
    end

    function Base.getindex(p::ParticlePool, tag::Symbol, idxs::Vector{Int})
        return p.data[tag][idxs]
    end

    function Base.getindex(p::ParticlePool, tag::Symbol, idx::Int)
        return p.data[tag][idx]
    end

    function Base.getindex(p::ParticlePool, re::Array{Bool,1})
        return p.pidx[re]
    end

    function Base.setindex!(p::ParticlePool, values::Array{Float64,1}, tag::Symbol)
        p.data[tag][p.pidx] = values
    end

    function Base.setindex!(p::ParticlePool, values::Array{Float64,1}, tag::Symbol, idxs::Array{Int64,1})
        p.data[tag][idxs] = values
    end

    function Base.setindex!(p::ParticlePool, value::Float64, tag::Symbol, idx::Int64)
        p.data[tag][idx] = value
    end

    function Base.display(p::ParticlePool)
        println("   $(p.tag),  ID: $(p.id),  nop: $(p.num_of_particles),  noa: $(length(p.links))")
    end

    function Base.keys(p::ParticlePool)
        return keys(p.data)
    end

    #------------------------------------------------------------------------------
    # set variables
    function set_var(p::ParticlePool, var_tags::Vector{Symbol})
        @assert length(keys(p)) == 0 "ERROR: Datatags Exist Already"
        foreach(var_tag->set_var(p, var_tag, 0.0), var_tags)
    end

    function set_var(p::ParticlePool, var_tag::Symbol, var_data::Vector{Float64})
        @assert isequal(length(var_data), p.num_of_particles) "ERROR: Vec Size != NOP"
        if var_tag in keys(p)
            p[var_tag] = var_data
            """for idx in 1:length(p.pidx)
            p[var_tag, p.pidx[idx]] = var_data[idx]
            end"""
        else
            get!(p.data, var_tag, var_data)
            append!(p.data[var_tag], zeros(p.num_of_reservation))
        end
    end

    function set_var(p::ParticlePool, var_tag::Symbol, var_data::Float64)
        set_var(p, var_tag, ones(p.num_of_particles).*var_data)
    end

    # Soft Remove and add particles in pool(within)
    function add(p::ParticlePool, indice::Vector{Int})
        # Soft add
        p.num_of_particles += length(indice)
        p.num_of_reservation -= length(indice)
        union!(p.pidx, indice)
        filter!(x->x∉indice, p.ridx)
        @assert isequal(length(p.pidx), p.num_of_particles) "ERROR: ADD p"
        @assert isequal(length(p.ridx), p.num_of_reservation) "ERROR: ADD r"
    end

    function remove(p::ParticlePool, indice::Vector{Int})
        p.num_of_particles -= length(indice)
        p.num_of_reservation += length(indice)
        union!(p.ridx, indice)
        filter!(x->x∉indice, p.pidx)
        @assert isequal(length(p.pidx), p.num_of_particles) "ERROR: REMOVE p"
        @assert isequal(length(p.ridx), p.num_of_reservation) "ERROR: REMOVE r"
    end

    # not enough reservation for new particles, memory re-allocating
    # DEBUG : optimize the append and merge scheme in future!
    function append(p::ParticlePool, size::Int)
        total_prev = p.num_of_particles+p.num_of_reservation
        append!(p.ridx, Array((total_prev+1):total_prev+size))
        for tag in keys(p)
            append!(p.data[tag], zeros(size))
        end
        p.num_of_reservation += size
    end

    # copy particles from src to dst pool and keep the particles in src
    function copy(dst::ParticlePool, src::ParticlePool, src_idx::Vector{Int})
        @assert isequal(length(keys(src)), length(keys(dst))) "ERROR: src NOV != dst NOV"
        if length(src_idx) > dst.num_of_reservation
            append(dst, length(src_idx)-dst.num_of_reservation)
        end
        for tag in keys(src)
            dst[tag, dst.ridx[1:length(src_idx)]] = src[tag, src_idx]
        end
"""        for idx in 1:length(src_idx)
            for tag in keys(src)
                dst[tag, dst.ridx[idx]] = src[tag, src_idx[idx]]
            end
        end"""
        add(dst, dst.ridx[1:length(src_idx)])
    end

    function copy(dst::ParticlePool, src::ParticlePool)
        copy(dst, src, src.pidx)
    end

    # copy particles from src to dst pool and remove the particles in src
    function transfer(dst::ParticlePool, src::ParticlePool, src_idx::Vector{Int})
        @assert !isequal(src, dst) "ERROR: Destination Pool == Source Pool"
        copy(dst, src, src_idx)
        remove(src, src_idx)
    end

    # merge src pool into dst
    function merge(dst::ParticlePool, src::ParticlePool)
        transfer(dst, src, src.pidx)
    end

    # soft clear a pool
    function clear(p::ParticlePool)
        if p.num_of_particles > 0
            remove(p, p.pidx)
        end
    end

    # clear a pool
    function empty(p::ParticlePool)
        empty!(p.pidx)
        empty!(p.ridx)
        p.num_of_particles = 0
        p.num_of_reservation = 0
    end

    # -----------------------------------------------------------------------------
    ## pool Boundary Functions
    function link(p1::ParticlePool, p2::ParticlePool, bd::GeoLib.Boundary)
        push!(p1.links, p2)
        push!(p2.links, p1)
        push!(p1.boundarys, GeoLib.opposite(bd))
        push!(p2.boundarys, bd)
    end

    function link(ps::Vector{Vector{ParticlePool}}, bd::Vector{GeoLib.Boundary})
        @assert ((length(ps)-1)==length(bd)) "ERROR: Boundary Num"
        for idx in 1:length(bd)
            link(ps[idx], ps[idx+1], bd[idx])
        end
    end

    function link(pg1::Vector{ParticlePool}, pg2::Vector{ParticlePool}, bd::GeoLib.Boundary)
        if length(pg1)==length(pg2)
            for idx in 1:length(pg1)
                link(pg1[idx], pg2[idx], bd)
            end
        elseif length(pg1)==1
            for idx in 1:length(pg2)
                link(pg1[1], pg2[idx], bd)
            end
        elseif length(pg2)==1
            for idx in range(1, length(pg1))
                link(pg1[idx], pg2[1], bd)
            end
        else
            println("ERROR: pool Number Not Match!")
        end
    end

    # ------------------------------------------------------------------------------
    #DEBUG Only allow same type pool and Inlet/Outlet transfer
    function transfer(dst::ParticlePool, src::ParticlePool, boundary::GeoLib.Boundary)
        # @assert src.tag == :Fluid, only Fluid pool excutes Transfer
        re::Array{Bool,1} = fill(false, src.num_of_particles)
        GeoLib.boundary_judge(boundary, src[:x], src[:y], src[:z], re)
        if dst.tag == src.tag
            transfer(dst, src, src[re])
        else
            if dst.tag == :Inlet
                println("Inlet leak: $(length(src[re]));  Remaining: $(length(src.pidx)-length(src[re]))")
                remove(src, src[re])
            elseif dst.tag == :Outlet
                println("Outlet leak: $(length(src[re]));  Remaining: $(length(src.pidx)-length(src[re]))")
                remove(src, src[re])
            elseif dst.tag == :Solid
                println("Wall $(dst.id) Penetration|Reflection: $(length(src[re]));  Remaining: $(length(src.pidx)-length(src[re]))")
                reflect(src, src[re], boundary)
            end
        end
    end

    function reflect(src::SPHParticles.ParticlePool, idx::Vector{Int}, boundary::GeoLib.Boundary)
"""        for idx in src.pidx
            if src[:x, idx] < 0.0
                src[:x, idx] = -src[:x, idx]
                src[:vx, idx] = -src[:vx, idx]
            elseif src[:x, idx] > 1.0
                src[:x, idx] = 2.0-src[:x, idx]
                src[:vx, idx] = -src[:vx, idx]
            end
            if src[:y, idx] < 0.0
                src[:y, idx] = -src[:y, idx]
                src[:vy, idx] = -src[:vy, idx]
            elseif src[:y, idx] > 1.0
                src[:y, idx] = 2.0-src[:y, idx]
                src[:vy, idx] = -src[:vy, idx]
            end
        end"""
        if length(idx)==0 return end
        #remove(src,idx)
        #println(typeof(boundary), length(idx))
        coeff = 1.0
        #println(typeof(boundary), " ",typeof(boundary) == GeoLib.PlateX," ",typeof(boundary) == GeoLib.PlateY)
        if typeof(boundary) == GeoLib.Boundary{GeoLib.PlateX}
            #println("before: ",src[:x, idx]," ",src[:y, idx]," ",src[:vx, idx]," ",src[:vy, idx])
            src[:vx, idx] = -src[:vx, idx]*coeff
            src[:vy, idx] = src[:vy, idx]*coeff
            src[:x, idx] = -src[:x, idx].+2.0*boundary.bound.cordin
            #println("after : ",src[:x, idx]," ",src[:vx, idx]," ",src[:vy, idx])
        elseif typeof(boundary) == GeoLib.Boundary{GeoLib.PlateY}
            #println("before: ",src[:x, idx]," ",src[:y, idx]," ",src[:vx, idx]," ",src[:vy, idx])
            src[:vx, idx] = src[:vx, idx]*coeff
            src[:vy, idx] = -src[:vy, idx]*coeff
            src[:y, idx] = -src[:y, idx].+2.0*boundary.bound.cordin
            #println("after : ",src[:x, idx]," ",src[:y, idx]," ",src[:vx, idx]," ",src[:vy, idx])
        end
    end

    function update(p::ParticlePool)
        if p.tag == :Solid return end
        if p.tag == :Outlet return end
        if p.tag == :Inlet
            #foreach(link->println(link.tag,"  ",link.id," ",length(link.ridx)," ",length(p.pidx)),p.links)
            foreach(dst->copy(dst, p), p.links)
        else
            #println(p.tag,"  ",length(p.links))
            #foreach(link->println(link.tag,"  ",link.id),p.links)
            foreach(idx->transfer(p.links[idx], p, p.boundarys[idx]), 1:length(p.links))
        end
    end

    ## Validation and Info Functions
    function check(p::ParticlePool)
        @assert length(p.links)==length(p.boundarys) " ERROR: pool boundary num: $(p.id)!"
    end
