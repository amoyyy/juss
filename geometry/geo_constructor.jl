    # position generating
    function rand(range::Vector{Float64}, num::Int)
        @assert range[2]>=range[1] "ERROR: rand position generating!"
        upper = range[2]
        lower = range[1]
        return Base.rand(num).*(upper-lower) .+ lower
    end

    function regular(xrange::Vector{Float64}, yrange::Vector{Float64}, ndims::Vector{Int})
        @assert xrange[2]>=xrange[1] && yrange[2]>=yrange[1]  "ERROR: regular position generating!"
        xarray = Array(range(xrange[1], length=ndims[1], stop=xrange[2]))
        yarray = Array(range(yrange[1], length=ndims[2], stop=yrange[2]))
        return meshgrid(xarray, yarray, ndims)
    end

    function regular(xrange::Vector{Float64}, yrange::Vector{Float64}, zrange::Vector{Float64}, ndims::Vector{Int})
        @assert xrange[2]>=xrange[1] && yrange[2]>=yrange[1] && zrange[2]>=zrange[1]  "ERROR: regular position generating!"
        xarray = Array(range(xrange[1], length=ndims[1], stop=xrange[2]))
        yarray = Array(range(yrange[1], length=ndims[2], stop=yrange[2]))
        zarray = Array(range(zrange[1], length=ndims[3], stop=zrange[2]))
        return meshgrid(xarray, yarray, zarray, ndims)
    end

    function rand_in_circle(center_point::Vector{Float64}, rrange::Vector{Float64}, num::Int)
        @assert rrange[2]>=rrange[1] "ERROR: rand circle position generating!"
        r::Vector{Float64} = rand(rrange, num)
        theta::Vector{Float64} = rand([0.0, 360.0], num)
        x::Vector{Float64} = r.*cosd.(theta).+center_point[1]
        y::Vector{Float64} = r.*sind.(theta).+center_point[2]
        return x, y
    end

    function regular_in_circle(center_point::Vector{Float64}, rrange::Vector{Float64}, ndims::Vector{Int})
        @assert rrange[2]>=rrange[1]  "ERROR: regular circle position generating!"
        r::Vector{Float64} = Array(range(rrange[1], length=ndims[1], stop=rrange[2]))
        theta::Vector{Float64} = Array(range(0.0, length=ndims[2], stop=360.0))
        tmp = meshgrid(r, theta, ndims)
        x::Vector{Float64} = tmp[1].*cosd.(tmp[2]).+center_point[1]
        y::Vector{Float64} = tmp[1].*sind.(tmp[2]).+center_point[2]
        return x, y
    end

    function rand_in_cylinder(center_point::Vector{Float64}, rrange::Vector{Float64}, zrange::Vector{Float64}, num::Int)
        @assert rrange[2]>=rrange[1] && zrange[2]>=zrange[1] "ERROR: rand cylinder position generating!"
        r::Vector{Float64} = rand(rrange, num)
        theta::Vector{Float64} = rand([0.0, 360.0], num)
        z::Vector{Float64} = rand(zrange, num)
        x::Vector{Float64} = r.*cosd.(theta).+center_point[1]
        y::Vector{Float64} = r.*sind.(theta).+center_point[2]
        return x, y, z
    end

    function regular_in_cylinder(center_point::Vector{Float64}, rrange::Vector{Float64}, zrange::Vector{Float64}, ndims::Vector{Int})
        @assert rrange[2]>=rrange[1] && zrange[2]>=zrange[1] "ERROR: regular cylinder position generating!"
        tmp_xy = regular_in_circle(center_point, rrange, ndims[1:2])
        x::Vector{Float64} = multi_array(tmp_xy[1], ndims[3], "Full")
        y::Vector{Float64} = multi_array(tmp_xy[2], ndims[3], "Full")
        z::Vector{Float64} = multi_array(Array(range(zrange[1],length=ndims[3],stop=zrange[2])), ndims[1]*ndims[2],"Single")
        return x, y, z
    end

    function rand_in_tank(xrange::Vector{Float64}, yrange::Vector{Float64}, zrange::Vector{Float64}, num::Int)
        @assert xrange[2]>=xrange[1] && yrange[2]>=yrange[1] && zrange[2]>=zrange[1] "ERROR: rand tank position generating!"
        x::Vector{Float64} = rand(xrange, num)
        y::Vector{Float64} = rand(yrange, num)
        z::Vector{Float64} = rand(zrange, num)
        return x, y, z
    end

    # support 2d mesh
    function meshgrid(xarray::Vector{Float64}, yarray::Vector{Float64}, ndims::Vector{Int})
        x::Vector{Float64} = zeros(ndims[1]*ndims[2])
        y::Vector{Float64} = zeros(ndims[1]*ndims[2])
        for i in 1:ndims[1]
            x[ndims[2]*(i-1)+1:ndims[2]*i] .= xarray[i]
            y[ndims[2]*(i-1)+1:ndims[2]*i] = yarray
        end
        return x, y
    end

    # support 3d mesh
    function meshgrid(xarray::Vector{Float64}, yarray::Vector{Float64}, zarray::Vector{Float64}, ndims::Vector{Int})
        x::Vector{Float64} = zeros(ndims[1]*ndims[2]*ndims[3])
        y::Vector{Float64} = zeros(ndims[1]*ndims[2]*ndims[3])
        z::Vector{Float64} = zeros(ndims[1]*ndims[2]*ndims[3])
        for i in 1:ndims[1]
            x[(i-1)*ndims[2]*ndims[3]+1:i*ndims[2]*ndims[3]] .= xarray[i]
            for j in 1:ndims[2]
                y[(i-1)*ndims[2]*ndims[3]+(j-1)*ndims[3]+1:(i-1)*ndims[2]*ndims[3]+j*ndims[3]] .= yarray[j]
                z[(i-1)*ndims[2]*ndims[3]+(j-1)*ndims[3]+1:(i-1)*ndims[2]*ndims[3]+j*ndims[3]] = zarray
            end
        end
        return x, y, z
    end

    # Mainly used in 3d mode
    function multi_array(base::Vector{Float64}, multi::Int, multi_mode::String)
        ret = zeros(length(base)*multi)
        if isequal(multi_mode, "Full")
            for i in Vector(1:multi)
                ret[Vector((i-1)*length(base)+1:i*length(base))] = base
            end
        elseif isequal(multi_mode, "Single")
            for ii in Vector(1:length(base))
                for jj in Vector(1:multi)
                    ret[(ii-1)*multi+1:ii*multi] = ones(multi).*base[ii]
                end
            end
        end
        return ret
    end
