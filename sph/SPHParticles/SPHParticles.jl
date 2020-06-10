###############################################################################
## SPHParticles
###############################################################################
module SPHParticles
    using NetCDF
    using GeoLib

    include("particle_pool.jl")
    include("nnp.jl")
    include("inlet_outlet.jl")
    include("particle_pool_union.jl")
    include("output_to_file.jl")
end
###############################################################################
