using Base.Threads
using NetCDF
using HDF5
using PyPlot
################################################################################
## Functions
################################################################################
function_list = ["clear_all(dir::String)"
                 "clear(dir::String, tag::String=\".png\")"
                 "gettime(dir::String)"
                 "display(dir::String, tag::String=\"Fluid\")"
                 "plotxy(dir::String)"
                 "plotxy_hdf5(dir::String)"
                 "plotxyz(dir::String)"
                 "plotT(dir::String, tmin::Real, tmax::Real, Tmin::Real, Tmax::Real)"]

for fuc in function_list
    println(fuc)
end
################################################################################
## 1
################################################################################
    function clear_all(dir::String)
        if isdir(dir)
            rm(dir, recursive=true)
            println("Directory: \"$dir\" removed!")
        else
            println("Directory: \"$dir\" doesn't exist!")
        end
    end

    function clear(dir::String, tag::String=".png")
        if isdir(dir)
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    if (occursin(tag, file))
                        filename = joinpath(root, file)
                        rm(filename)
                    end
                end
            end
            println("$tag Files in Directory: \"$dir\" removed!")
        else
            println("Directory: \"$dir\" doesn't exist!")
        end
    end

################################################################################
## 2
################################################################################
    function gettime(dir::String)
        time_vec::Vector{Float64} = Float64[]
        if isdir(dir)
            for (root, dirs, files) in walkdir(dir)
                for dir in dirs
                    push!(time_vec, parse(Float64,(dir[3:end-1])))
                end
            end
        end
        return sort!(time_vec)
    end

    function getstep(dir::String)
        step_vec::Vector{Float64} = Float64[]
        if isdir(dir)
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    if (occursin(".hdf5", file))
                        istart = findfirst("_", file)[1]
                        iend = findfirst(".hdf5", file)[1]
                        push!(step_vec, parse(Int64,(file[istart+1:iend-1])))
                    end
                end
            end
        end
        return sort!(step_vec)
    end

    # Support Function
    function get_particles_num_onetime(timedir::String)
        for (root, dirs, files) in walkdir(timedir)
            for file in files
                if (occursin(".nc", file))
                    filename = joinpath(root, file)
                    return length(ncread(filename, "x"))
                end
            end
        end
    end

    # Support Function
    function get_value(timedir::String, var_name::String, tag::String="")
        num = get_particles_num_onetime(timedir)
        var_all::Vector{Float64} = Float64[]
        for (root, dirs, files) in walkdir(timedir)
            for file in files
                if occursin(".nc", file) && occursin(tag, file)
                    filename = joinpath(root, file)
                    var::Vector{Float64} = ncread(filename, var_name)
                    var_all = [var_all; var]
                    ncclose(filename)
                end
            end
        end
        return var_all
    end

    function display_data(dir::String, tag::String="Fluid", step::Int64=1)
        times = gettime(dir)
        times = times[1:step:length(times)]
        io = IOBuffer()
        for str in  ["time", "num_of_particles", "sum(m)", "sum(V)", "sum(px)", "sum(py)", "sum(pz)", "sum(Ev)", "sum(m_re)","sum(px_re)", "sum(py_re)", "sum(pz_re)",
                    "max(rho)", "min(rho)", "max(P)", "min(P)", "max(x)", "min(x)", "max(y)", "min(y)","max(z)", "min(z)", "max(ax)", "max(ay)","max(az)","max(T)", "min(T)", "ave(T)"]
            print(io,str,",")
        end
        println(io)
        for time in times
            println("reading data: time=$time")
            timedir = joinpath(dir, string("t=$time","s"))
            m = get_value(timedir, "m", tag)
            rho = get_value(timedir, "rho", tag)
            V = 1.0./get_value(timedir, "V", tag)
            P = get_value(timedir, "p", tag)
            x = get_value(timedir, "x", tag)
            y = get_value(timedir, "y", tag)
            z = get_value(timedir, "z", tag)
            vx = get_value(timedir, "vx", tag)
            vy = get_value(timedir, "vy", tag)
            vz = get_value(timedir, "vz", tag)
            ax = get_value(timedir, "ax", tag)
            ay = get_value(timedir, "ay", tag)
            az = get_value(timedir, "az", tag)
            T = get_value(timedir, "T", tag)

            px = m.*vx
            py = m.*vy
            pz = m.*vz
            ev = 0.5*m.*(vx.*vx+vy.*vy+vz.*vz)
            m_re = V.*rho
            px_re = m_re.*vx
            py_re = m_re.*vy
            pz_re = m_re.*vz

            println(io,time,",",length(m),",",
                    sum(m),",", sum(V),",", sum(px),",",sum(py),",",sum(pz),",",sum(ev),",",
                    sum(m_re),",", sum(px_re),",", sum(py_re),",",sum(pz_re),",",
                    maximum(rho),",", minimum(rho),",", maximum(P),",", minimum(P),",",
                    maximum(x),",", minimum(x),",", maximum(y),",", minimum(y),",",maximum(z),",", minimum(z),",",
                    maximum(ax),",", maximum(ay),",",maximum(az),",",
                    maximum(T),",", minimum(T),",", sum(T)/length(T)
                    )
        end
        OUTPUT_FILE = joinpath(dir, "Values_$tag.csv")
        open(OUTPUT_FILE, "w") do csvio
            buffer = String(take!(io))
            write(csvio, buffer)
        end
        println("Output File: $OUTPUT_FILE ")
    end

    function display_t(dir::String, tag::String="Fluid", step::Int64=1)
        header_print::Bool = true
        times = gettime(dir)
        times = times[1:step:length(times)]
        io = IOBuffer()
        for time in times
            println("reading data: time=$time")
            timedir = joinpath(dir, string("t=$time","s"))
            x = get_value(timedir, "x", tag)
            y = get_value(timedir, "y", tag)
            z = get_value(timedir, "z", tag)
            T = get_value(timedir, "T", tag)
            y_level = minimum(y[y.>0.5*(minimum(y)+maximum(y))])
            idx = Array(1:length(x))
            idx = idx[y.==y_level]
            #x_new = x[idx]
            #t_new = t[idx]
            if header_print
                print(io,"y=$y_level|x=,")
                for i in idx
                    print(io,x[i],",")
                end
                println(io,"")
                header_print = false
            end
            print(io,time,",")
            for i in idx
                print(io,T[i],",")
            end
            println(io,"")
        end
        OUTPUT_FILE = joinpath(dir, "Temperature_Distribution_$tag.csv")
        open(OUTPUT_FILE, "w") do csvio
            buffer = String(take!(io))
            write(csvio, buffer)
        end
        println("Output File: $OUTPUT_FILE ")
    end

################################################################################
## 3
################################################################################
    function plot(filename::String)
        if (occursin(".nc", filename) && isfile(filename))
            println(filename)
            x::Array{Float64,1} = ncread(filename, "x")
            y::Array{Float64,1} = ncread(filename, "y")
            z::Array{Float64,1} = ncread(filename, "z")
            ncclose(filename)
            PyPlot.scatter(x,y,z)
            PyPlot.savefig(replace(filename, Pair(".nc", ".png")))
            PyPlot.clf()
        else
            println("$filename doesn't exist!")
        end
    end

    function plotxy(dir::String, xlims::Vector{Float64}, ylims::Vector{Float64})
        times = gettime(dir)
        PyPlot.figure()
        #@threads
        for time in times
            println("plotting: time=$time")
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    filename = joinpath(root, file)
                    if occursin(".nc", filename) && occursin(string("t=",time,"s"), filename)
                        x::Vector{Float64} = ncread(filename, "x")
                        y::Vector{Float64} = ncread(filename, "y")
                        if occursin("Solid", file)
                            PyPlot.plot(x,y,"k.")
                        elseif occursin("Fluid", file)
                            PyPlot.plot(x,y,"b.")
                        end
                        ncclose(filename)
                    end
                end
            end
            PyPlot.xlim(xlims[1],xlims[2])
            PyPlot.ylim(ylims[1],ylims[2])
            PyPlot.title("t = $time s")
            PyPlot.savefig(string(joinpath(dir,"t=$time"), "s.png"), dpi=600)
            PyPlot.clf()
        end
        PyPlot.close()
    end

    function plotxy_hdf5(dir::String, xlims::Vector{Float64}, ylims::Vector{Float64})
        steps = getstep(dir)
        PyPlot.figure()
        for step in steps
            println("plotting: step=$step")
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    filename = joinpath(root, file)
                    if occursin(string("_", Int64(step),".hdf5"), filename)
                        fid = h5open(filename, "r")
                        data = read(fid, "particles")
                        x::Vector{Float64} = data["fluid"]["arrays"]["x"]
                        y::Vector{Float64} = data["fluid"]["arrays"]["y"]
                        x1::Vector{Float64} = data["solid"]["arrays"]["x"]
                        y1::Vector{Float64} = data["solid"]["arrays"]["y"]
                        PyPlot.plot(x1,y1,"k.")
                        PyPlot.plot(x,y,"b.")
                        close(fid)
                    end
                end
            end
            PyPlot.xlim(xlims[1],xlims[2])
            PyPlot.ylim(ylims[1],ylims[2])
            PyPlot.title("step = $step")
            PyPlot.savefig(string(joinpath(dir,"step=$step"), ".png"), dpi=600)
            PyPlot.clf()
        end
        PyPlot.close()
    end

    function plotxyz(dir::String)
        times = gettime(dir)
        PyPlot.figure()
        for time in times
            println("plotting: time=$time")
            for (root, dirs, files) in walkdir(dir)
                for file in files
                    filename = joinpath(root, file)
                    if occursin(".nc", filename) && occursin(string("t=",time,"s"), filename)
                        x::Vector{Float64} = ncread(filename, "x")
                        y::Vector{Float64} = ncread(filename, "y")
                        z::Vector{Float64} = ncread(filename, "z")
                        if occursin("Solid", file)
                            PyPlot.scatter3D(x,y,z,s=10,"k")
                        elseif occursin("Fluid", file)
                            PyPlot.scatter3D(x,y,z,s=10,"blue")
                        elseif occursin("Inlet", file)
                            PyPlot.scatter3D(x,y,z,s=10,"y")
                        end
                        ncclose(filename)
                    end
                end
            end
            PyPlot.title("t = $time s")
            PyPlot.savefig(string(joinpath(dir,"t=$time"), "s.png"), dpi=600)
            PyPlot.clf()
        end
        PyPlot.close()
    end

    function plotT(dir::String, tmin::Real, tmax::Real, Tmin::Real, Tmax::Real)
        times = gettime(dir)
        cm = PyPlot.get_cmap("rainbow")
        PyPlot.figure()
        for time in times
            if tmin<time<tmax
                println("plotting: time=$time")
                x::Vector{Float64} = Float64[]
                y::Vector{Float64} = Float64[]
                z::Vector{Float64} = Float64[]
                t::Vector{Float64} = Float64[]
                for (root, dirs, files) in walkdir(dir)
                    for file in files
                        filename = joinpath(root, file)
                        if occursin(".nc", filename) && occursin(string("t=",time,"s"), filename)
                            x = [x;ncread(filename, "x")]
                            y = [y;ncread(filename, "y")]
                            z = [z;ncread(filename, "z")]
                            t = [t;ncread(filename, "T")]
                            ncclose(filename)
                        end
                    end
                end
                fig = PyPlot.scatter(x, y, c=t, s=10, vmin=Tmin, vmax=Tmax, cmap=cm, label=string(time))
                PyPlot.colorbar(fig)
                PyPlot.title("t = $time s")
                PyPlot.savefig(string(joinpath(dir,"t=$time"), "s.png"), dpi=600)
                PyPlot.clf()
            end
        end
        PyPlot.close()
    end
