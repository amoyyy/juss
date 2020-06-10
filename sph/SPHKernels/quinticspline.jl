    struct QuinticSpline <: SPHKernels.Kernel
        dim::Int
        factor::Float64
        radius_scale::Float64
        QuinticSpline(dim, factor, radius_scale) = new(dim, factor, radius_scale)
        QuinticSpline(dim) =
        QuinticSpline(dim, map(x->begin
                            if x == 1
                                factor = 1.0 / 120.0
                            elseif x == 2
                                factor = 7.0 / 478.0 / pi
                            elseif x == 3
                                factor = 3.0 / 359.0 / pi
                            end
                            return factor
                        end,dim), 3.0)
    end

    function get_coef(p::SPHKernels.QuinticSpline)
        # The inflection points for the polynomial are obtained as
        # http://www.wolframalpha.com/input/?i=%28%283-x%29%5E5+-+6*%282-x%29%5E5+%2B+15*%281-x%29%5E5%29%27%27
        # the only permissible value is taken
        return 0.759298480738450
    end

    function kernel(p::SPHKernels.QuinticSpline, rij::Float64=1.0, h::Float64=1.0)
        h1 = 1.0 / h
        q = rij * h1

        # get the kernel normalizing factor
        if p.dim == 1
            fac = p.factor * h1
        elseif p.dim == 2
            fac = p.factor * h1 * h1
        elseif p.dim == 3
            fac = p.factor * h1 * h1 * h1
        end

        tmp3 = 3.0 - q
        tmp2 = 2.0 - q
        tmp1 = 1.0 - q

        val = 0.0
        if (q > 3.0)
            val = 0.0
        elseif (q > 2.0)
            val = tmp3 * tmp3 * tmp3 * tmp3 * tmp3
        elseif (q > 1.0)
            val = tmp3 * tmp3 * tmp3 * tmp3 * tmp3
            val -= 6.0 * tmp2 * tmp2 * tmp2 * tmp2 * tmp2
        else
            val = tmp3 * tmp3 * tmp3 * tmp3 * tmp3
            val -= 6.0 * tmp2 * tmp2 * tmp2 * tmp2 * tmp2
            val += 15.0 * tmp1 * tmp1 * tmp1 * tmp1 * tmp1
        end

        return val * fac
    end

    function dwdq(p::SPHKernels.QuinticSpline, rij::Float64=1.0, h::Float64=1.0)
        h1 = 1.0 / h
        q = rij * h1

        # get the kernel normalizing factor
        if p.dim == 1
            fac = p.factor * h1
        elseif p.dim == 2
            fac = p.factor * h1 * h1
        elseif p.dim == 3
            fac = p.factor * h1 * h1 * h1
        end

        tmp3 = 3.0 - q
        tmp2 = 2.0 - q
        tmp1 = 1.0 - q

        val = 0.0
        # compute the gradient
        if (rij > 1e-12)
            if (q > 3.0)
                val = 0.0
            elseif (q > 2.0)
                val = -5.0 * tmp3 * tmp3 * tmp3 * tmp3
            elseif (q > 1.0)
                val = -5.0 * tmp3 * tmp3 * tmp3 * tmp3
                val += 30.0 * tmp2 * tmp2 * tmp2 * tmp2
            else
                val = -5.0 * tmp3 * tmp3 * tmp3 * tmp3
                val += 30.0 * tmp2 * tmp2 * tmp2 * tmp2
                val -= 75.0 * tmp1 * tmp1 * tmp1 * tmp1
            end
        end
        return val * fac
    end

    function gradient(p::SPHKernels.QuinticSpline, rij=1.0, h=1.0, xij=[0., 0., 0.], grad=[0., 0., 0.])
        h1 = 1. / h
        tmp = 0.0
        wdash = 0.0

        # compute the gradient.
        if (rij > 1e-12)
            wdash = dwdq(p, rij, h)
            tmp = wdash * h1 / rij
        end

        for i in 1:p.dim
            grad[i] = tmp * xij[i]
        end
    end

    function gradient_h(p::SPHKernels.QuinticSpline, rij::Float64=1.0, h::Float64=1.0, xij=[0., 0., 0.])
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

        tmp3 = 3.0 - q
        tmp2 = 2.0 - q
        tmp1 = 1.0 - q

        # compute the kernel & gradient at q
        if (q > 3.0)
            w = 0.0
            dw = 0.0
        elseif (q > 2.0)
            w = tmp3 * tmp3 * tmp3 * tmp3 * tmp3
            dw = -5.0 * tmp3 * tmp3 * tmp3 * tmp3
        elseif (q > 1.0)
            w = tmp3 * tmp3 * tmp3 * tmp3 * tmp3
            w -= 6.0 * tmp2 * tmp2 * tmp2 * tmp2 * tmp2
            dw = -5.0 * tmp3 * tmp3 * tmp3 * tmp3
            dw += 30.0 * tmp2 * tmp2 * tmp2 * tmp2
        else
            w = tmp3 * tmp3 * tmp3 * tmp3 * tmp3
            w -= 6.0 * tmp2 * tmp2 * tmp2 * tmp2 * tmp2
            w += 15. * tmp1 * tmp1 * tmp1 * tmp1 * tmp1
            dw = -5.0 * tmp3 * tmp3 * tmp3 * tmp3
            dw += 30.0 * tmp2 * tmp2 * tmp2 * tmp2
            dw -= 75.0 * tmp1 * tmp1 * tmp1 * tmp1
        end

        return -fac * h1 * (dw * q + w * p.dim)
    end
