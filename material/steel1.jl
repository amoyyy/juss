################################################################################
    struct CFRIHXSteel <: Material end
    #密度(kg/m3)
    #ro = (mat::Steel, tc::Float64, pressure::Float64) -> 7989 + (0.127 + 1.51e-5*tc*tc)*tc
    function ro(mat::CFRIHXSteel, tc::Float64, pressure::Float64)
        return 8000.0
    end

    #比热(J/kg*K)
    #c_p = (mat::Steel, tc::Float64, pressure::Float64) -> 500 + (0.072 - (6.37e-4 + 1.73e-6*tc)*tc)*tc
    function c_p(mat::CFRIHXSteel, tc::Float64, pressure::Float64)
        #return 533.27
        return 575.0
    end

    #导热系数(W/m*K)
    #lamda = (mat::Steel, tc::Float64, pressure::Float64) -> 10.679 + 1.18e-2*tc
    function lamda(mat::CFRIHXSteel, tc::Float64, pressure::Float64)
        return 20.25
    end
