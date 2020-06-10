    function f_step(in::Float64)
        if in<0.0
            return 0.0
        else
            return 1.0
        end
    end

    function f_delta(in::Float64)
        if in==0.0
            return 1.0
        else
            return 0.0
        end
    end
