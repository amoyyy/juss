    # The output action is now excuted in SPHSolver.solve, the artitechture is to changed later
    function output(p::SPHParticles.ParticlePool, dirpath::String, filename::String)
        indexatts = Dict("longname" => "id",
                        "units"    => "1")
        xatts = Dict("longname" => "X",
                       "units"    => "m")
        yatts = Dict("longname" => "Y",
                        "units"    => "m")
        zatts = Dict("longname" => "Z",
                        "units"    => "m")
        timeatts = Dict("longname" => "time",
                        "units"    => "s")

        (!isdir(dirpath)) && mkpath(dirpath)
        filepath = joinpath(dirpath, filename)
        isfile(filepath) && rm(filepath)

        for tag in keys(p)
            varname = String(tag)
            varatts = Dict("longname" => varname)
            nccreate(filepath, varname, "id", Array(1:p.num_of_particles), indexatts, atts=varatts)
            ncwrite(p[tag], filepath, varname)
        end
        ncclose(filepath)
    end
