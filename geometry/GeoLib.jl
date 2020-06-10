################################################################################
## Basic Geometry Boundary Type
################################################################################
module GeoLib
    abstract type AbstractBoundary end
    # Plate
    mutable struct PlateX <: AbstractBoundary
        cordin::Float64
        PlateX(cordin=0.0) = new(cordin)
    end

    mutable struct PlateY <: AbstractBoundary
        cordin::Float64
        PlateY(cordin=0.0) = new(cordin)
    end

    mutable struct PlateZ <: AbstractBoundary
        cordin::Float64
        PlateZ(cordin=0.0) = new(cordin)
    end

    ## Curve
    mutable struct CylinderZ <: AbstractBoundary
        center_cordin::Array{Float64}
        radius::Float64
        radius_square::Float64
        CylinderZ(center_cordin::Vector{Float64}=[0.0,0.0], radius::Float64=0.0) = new(center_cordin, radius, radius^2)
    end

    struct NoBoundary <: AbstractBoundary end

    function boundary_judge(go::AbstractBoundary, x::Array{Float64,1}, re::Array{Bool,1}, op::Function)
    end

    function boundary_judge(go::PlateX, xcor::Array{Float64,1}, ycor::Array{Float64,1}, zcor::Array{Float64,1}, re::Array{Bool,1}, op::Function)
        map!(x->(op(x,go.cordin)), re, xcor)
    end

    function boundary_judge(go::PlateY, xcor::Array{Float64,1}, ycor::Array{Float64,1}, zcor::Array{Float64,1}, re::Array{Bool,1}, op::Function)
        map!(y->(op(y, go.cordin)), re, ycor)
    end

    function boundary_judge(go::PlateZ, xcor::Array{Float64,1}, ycor::Array{Float64,1}, zcor::Array{Float64,1}, re::Array{Bool,1}, op::Function)
        map!(z->(op(z, go.cordin)), re, zcor)
    end

    function boundary_judge(go::CylinderZ, cordin::Array{Float64,2}, re::Array{Bool,1}, op::Function)
        idx = Vector(1:length(cordin)/2)
        map!(x->(op(sum(abs2.(cordin[1:2,x].-go.center_cordin)), radius_square)), re, idx)
    end

    function boundary_judge(go::CylinderZ, xcor::Array{Float64,1}, ycor::Array{Float64,1}, zcor::Array{Float64,1}, re::Array{Bool,1}, op::Function)
        idx = Vector(1:length(xcor))
        map!(x->(op(abs2(xcor[x]-go.center_cordin[1])+abs2(ycor[x]-go.center_cordin[2]), go.radius_square)), re, idx)
    end

    function boundary_judge(go::NoBoundary, xcor::Array{Float64,1}, ycor::Array{Float64,1}, zcor::Array{Float64,1}, re::Array{Bool,1}, op::Function)
        re .= false
    end

################################################################################
## Define Area Type
################################################################################
    mutable struct Boundary{T <: AbstractBoundary}
        bound::T
        op::Function
        Boundary(bd=PlateX(0.0), op=<) = new{typeof(bd)}(bd, op)
    end

    function boundary_judge(b::Boundary, xcor::Array{Float64,1}, ycor::Array{Float64,1}, zcor::Array{Float64,1}, re::Array{Bool,1})
        boundary_judge(b.bound, xcor, ycor, zcor, re, b.op)
    end

################################################################################
## Basic Geometry Operators
################################################################################
    function opposite(bd::Boundary)
        if isequal(bd.op,<) || isequal(bd.op,isless)
            op = >
        elseif isequal(bd.op,>)
            op = <
        elseif isequal(bd.op,==)
            op = ==
        end
        return Boundary(bd.bound, op)
    end
################################################################################
    include("geo_constructor.jl")

end
