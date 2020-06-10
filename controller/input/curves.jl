    function flowrate_curve_constant(t::Float64)
        return 1.0
    end

    function flowrate_curve_1(t::Float64)
        return f_step(t)-0.9*(t-50.0)*(f_step(t-50.0)-f_step(t-51.0))-(f_step(t-51.0)-f_step(t-100.0))*0.9-(f_step(t-100.0)-f_step(t-101.0))*(101.0-t)*0.9
    end

    function flowrate_curve_2(t::Float64)
        return f_step(t)-(f_step(t-151.0)*0.9+(t-150.0)*(f_step(t-150.0)-f_step(t-151.0))*0.9)
    end

    function flowrate_sind(t::Float64)
        return 0.6+0.4*sin(2.0*pi*t/10.0)
    end
