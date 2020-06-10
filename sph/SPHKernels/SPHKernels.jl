module SPHKernels
    abstract type Kernel end

    include("cubicspline.jl")
    include("quinticspline.jl")
    include("gaussian.jl")
    include("supergaussian.jl")

    AllKernels = Union{Gaussian,SuperGaussian,CubicSpline,QuinticSpline}
end
