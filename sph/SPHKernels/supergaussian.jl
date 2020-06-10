
    struct SuperGaussian <: SPHKernels.Kernel
        dim::Int
        factor::Float64
        SuperGaussian(dim, factor) = new(dim, factor)
        SuperGaussian(dim) =
        SuperGaussian(dim, map(x->begin
                            factor = 1.0 / sqrt(pi)
                            if x>1
                                factor *= (1.0 / sqrt(pi))
                            end
                            if x>2
                                factor *= (1.0 / sqrt(pi))
                            end
                            return factor
                            end,dim))
    end

    function get_coef(p::SPHKernels.SuperGaussian)
        # Found inflection point using sympy.
        if p.dim == 1
            return 0.584540507426389
        elseif p.dim == 2
            return 0.6021141014644256
        else
            return 0.615369528365158
        end
    end

    function kernel(p::SPHKernels.SuperGaussian, rij::Float64=1.0, h::Float64=1.0)
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

        val = 0.0
        if (0.0 < q < 3.0)
            q2 = q * q
            val = exp(-q2) * (1.0 + p.dim * 0.5 - q2) * fac
        end

        return val
    end

    function dwdq(p::SPHKernels.SuperGaussian, rij::Float64=1.0, h::Float64=1.0)
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
        # compute the gradient
        val = 0.0
        if (q < 3.0)
            if (rij > 1e-12)
                q2 = q * q
                val = q * (2.0 * q2 - p.dim - 4) * exp(-q2)
            end
        end
        return val * fac
    end

    function gradient(p::SPHKernels.SuperGaussian, rij::Float64=1.0, h::Float64=1.0, xij=[0., 0., 0.], grad=[0, 0, 0])
        h1 = 1. / h
        # compute the gradient.
        if (rij > 1e-12)
            wdash = dwdq(prij, h)
            tmp = wdash * h1 / rij
        else
            tmp = 0.0
        end

        for i in 1:p.dim
            grad[i] = tmp * xij[i]
        end
    end

    function gradient_h(p::SPHKernels.SuperGaussian, rij::Float64=1.0, h::Float64=1.0, xij=[0., 0., 0.])
        h1 = 1. / h
        q = rij * h1
        d = p.dim

        # get the kernel normalizing factor
        if d == 1
            fac = p.factor * h1
        elseif d == 2
            fac = p.factor * h1 * h1
        elseif d == 3
            fac = p.factor * h1 * h1 * h1
        end

        # kernel and gradient evaluated at q
        val = 0.0
        if (q < 3.0)
            q2 = q * q
            val = (-d * d * 0.5 + 2.0 * d * q2 - d - 2.0 * q2 * q2 + 4 * q2) * exp(-q2)
        end

        return -fac * h1 * val
    end
