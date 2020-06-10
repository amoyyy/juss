

    # `CubicSpline`
    struct CubicSpline <: SPHKernels.Kernel
        dim::Int
        factor::Float64
        CubicSpline(dim, factor) = new(dim, factor)
        CubicSpline(dim) =
        CubicSpline(dim, map(x->begin
                            if x == 3
                                factor = 1.0 / pi
                            elseif x == 2
                                factor = 10.0 / 7.0 / pi
                            else
                                factor = 2.0 / 3.0
                            end
                            return factor
                            end,dim))

    end

    function get_coef(p::SPHKernels.CubicSpline)
        return 2. / 3
    end

    function kernel(p::SPHKernels.CubicSpline, rij::Float64=1.0, h::Float64=1.0)
        h1 = 1. / h
        q = rij * h1

        # get the kernel normalizing factor
        if p.dim == 1
            fac = p.factor * h1
        elseif p.dim == 2
            fac = p.factor * h1 * h1
        elseif p.dim == 3
            fac = p.factor * h1 * h1 * h1
        end

        tmp2 = 2.0 - q
        val = 0.0
        if (q > 2.0)
            val = 0.0
        elseif (q > 1.0)
            val = 0.25 * tmp2 * tmp2 * tmp2
        elseif (q > 0.0)
            val = 1 - 1.5 * q * q * (1.0 - 0.5 * q)
        end

        return val * fac
    end


    function dwdq(p::SPHKernels.CubicSpline, rij::Float64=1.0, h::Float64=1.0)
        h1 = 1. / h
        q = rij * h1

        # get the kernel normalizing factor ( sigma )
        if p.dim == 1
            fac = p.factor * h1
        elseif p.dim == 2
            fac = p.factor * h1 * h1
        elseif p.dim == 3
            fac = p.factor * h1 * h1 * h1
        end

        # compute sigma * dw_dq
        tmp2 = 2. - q
        if (rij > 1e-12)
            if (q > 2.0)
                val = 0.0
            elseif (q > 1.0)
                val = -0.75 * tmp2 * tmp2
            else
                val = -3.0 * q * (1 - 0.75 * q)
            end
        else
            val = 0.0
        end

        return val * fac
    end

    function gradient(p::SPHKernels.CubicSpline, rij::Float64=1.0, h::Float64=1.0, xij=[0., 0, 0], grad=[0, 0, 0])
        h1 = 1. / h
        # compute the gradient.
        if (rij > 1e-12)
            wdash = dwdq(p, rij, h)
            tmp = wdash * h1 / rij
        else
            tmp = 0.0
        end

        for i in 1:p.dim
            grad[i] = tmp * xij[i]
        end
    end

    function gradient_h(p::SPHKernels.CubicSpline, rij::Float64=1.0, h::Float64=1.0)
        h1 = 1. / h
        q = rij * h1

        # get the kernel normalizing factor
        if p.dim == 1
            fac = p.factor * h1
        elseif p.dim == 2
            fac = p.factor * h1 * h1
        elseif p.dim == 3
            fac = p.factor * h1 * h1 * h1
        end

        # kernel and gradient evaluated at q
        tmp2 = 2. - q
        if (q > 2.0)
            w = 0.0
            dw = 0.0
        elseif (q > 1.0)
            w = 0.25 * tmp2 * tmp2 * tmp2
            dw = -0.75 * tmp2 * tmp2
        else
            w = 1 - 1.5 * q * q * (1 - 0.5 * q)
            dw = -3.0 * q * (1 - 0.75 * q)
        end

        return -fac * h1 * (dw * q + w * p.dim)
    end
