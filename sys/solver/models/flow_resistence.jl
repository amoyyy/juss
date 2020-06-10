###############################################################################
## FlowResistence
###############################################################################
export s_flow_resistence_loss

"""
//For circular pipes;
//------------------------------------------------------------------------------
//  d:      hydraulic diameter
//  D:      pipe diameter(In circular pipes d=D)
//  area:   cross-sectional area

//  V:      velocity of the fluid
//  RHO:    density of the fulid
//  mu:
//  e:
//------------------------------------------------------------------------------
"""

    re = (d::Float64, V::Float64,  RHO::Float64,  eta::Float64) -> abs(V)*RHO*d/eta

    function s_flow_resistence_loss(V::Float64,  RHO::Float64,  d::Float64,  eta::Float64,  e::Float64)
        re = re(d, V, RHO, eta);
        f = s_friction_factor(re,e,d);
        return 0.5 * f * RHO * V * fabs(V) / d;    # N/m^3
    end

    function s_friction_factor(Re::Float64, e::Float64, d::Float64)
        fac::Float64 = 0
        if 64<=Re<2200
            fac = 64/Re
        elseif 2200<=Re<3000
            f1 = 0.029091;
            f2 = s_friction_factor(3000.0,e,d);
            weight = 3.75 - 8250.0 / Re;
            fac = (1-weight)*f1 + weight*f2;
        elseif Re>3000
            tmp1 = e/d;
            tmp2 = 1.14-2.0*log10(tmp1 + 21.25/(Re^0.9));
            tmp3 = -2.0*log10(tmp1/3.7+2.51/Re*tmp2);
            fac = 1.0 / (tmp3 * tmp3);
        elseif 0<=Re<64
            fac = NaN
        elseif Re<0
            nothing
        end
        return fac
    end
