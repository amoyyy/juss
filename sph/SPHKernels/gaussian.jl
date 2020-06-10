    """
    References
    ----------
    .. [Liu2010] `M. Liu, & G. Liu, Smoothed particle hydrodynamics (SPH)
        an overview and recent developments, "Archives of computational
        methods in engineering", 17.1 (2010), pp. 25-76.
        <http//link.springer.com/article/10.1007/s11831-010-9040-7>`_
    """
    struct Gaussian <: SPHKernels.Kernel
        dim::Int
        factor::Float64
        Gaussian(dim, factor) = new(dim, factor)
        Gaussian(dim) =
        Gaussian(dim, map(x->begin
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

    """
    References
    ----------
    .. [Monaghan1992] `J. Monaghan, Smoothed Particle Hydrodynamics, "Annual
        Review of Astronomy and Astrophysics", 30 (1992), pp. 543-574.
        <http://adsabs.harvard.edu/abs/1992ARA&A..30..543M>`_
    """

    ################################################################################
    # Functions
    ################################################################################
    function kernel(p::SPHKernels.Gaussian, rij::Float64=1.0, h::Float64=1.0)
        h1 = 1. / h
        q = rij * h1

        # get the kernel normalizing factor
        if p.dim == 1
            fac = p.factor * h1
        elseif p.dim == 2
            fac = p.factor * h1 * h1
        else p.dim == 3
            fac = p.factor * h1 * h1 * h1
        end

        val = 0.0
        if (q < 3.0)
            val = exp(-q * q) * fac
        end

        return val
    end

    function gradient(p::SPHKernels.Gaussian, rij::Float64=1.0, h::Float64=1.0, xij::Vector{Float64}=[0., 0., 0.], grad::Vector{Float64}=[0, 0, 0])
        h1 = 1. / h

        # compute the gradient.
        if (rij > 1e-12)
            wdash = dwdq(p,rij, h)
            tmp = wdash * h1 / rij
        else
            tmp = 0.0
        end

        for i in 1:p.dim
            grad[i] = tmp * xij[i]
        end
    end

    function gradient_h(p::SPHKernels.Gaussian, rij::Float64=1.0, h::Float64=1.0)
        h1 = 1. / h
        q = rij * h1

        # get the kernel normalizing factor
        if p.dim == 1
            fac = p.factor * h1
        elseif p.dim == 2
            fac = p.factor * h1 * h1
        else p.dim == 3
            fac = p.factor * h1 * h1 * h1
        end

        # kernel and gradient evaluated at q
        w = 0.0
        dw = 0.0
        if (q < 3.0)
            w = exp(-q * q)
            dw = -2.0 * q * w
        end
        return -fac * h1 * (dw * q + w * p.dim)
    end

    function dwdq(p::SPHKernels.Gaussian, rij::Float64=1.0, h::Float64=1.0)
        h1 = 1. / h
        q = rij * h1

        # get the kernel normalizing factor
        if p.dim == 1
            fac = p.factor * h1
        elseif p.dim == 2
            fac = p.factor * h1 * h1
        else p.dim == 3
            fac = p.factor * h1 * h1 * h1
        end

        # compute the gradient
        val = 0.0
        if (q < 3.0)
            if (rij > 1e-12)
                val = -2.0 * q * exp(-q * q)
            end
        end

        return val * fac
    end

    function get_coef(p::SPHKernels.Gaussian)
        return 0.70710678118654746
    end
