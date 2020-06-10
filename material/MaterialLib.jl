################################################################################
## MaterialLib
################################################################################
module MaterialLib
    abstract type Material end
    include("sodium.jl")
    include("water.jl")
    include("steel.jl")
    include("steel1.jl")

    function getMat(type::String)
        if type == "Na"
            return Sodium()
        elseif type == "H2O"
            return Water()
        elseif type == "Steel"
            return Steel()
        elseif type == "CFR-IHX-Steel"
            return CFRIHXSteel()
        end
    end

    function getMat(type::Symbol)
        return getMat(string(type))
    end

    #焓(J/kg)
    function enth(mat::Material, tc::Float64, pressure::Float64) end
    #密度(kg/m3)
    function ro(mat::Material, tc::Float64, pressure::Float64) end
    #比热(J/kg*K)
    function c_p(mat::Material, tc::Float64, pressure::Float64) end
    #动力粘度(kg/m*s)
    function eta(mat::Material, tc::Float64, pressure::Float64) end
    #导热系数(W/m*K)
    function lamda(mat::Material, tc::Float64, pressure::Float64) end
    # 材料声速
    function cs(mat::Material, tc::Float64, pressure::Float64) end
    # 反解温度
    function temperature(mat::Material, enthalpy::Float64, pressure::Float64) end
    #焓对压力偏导(J/kg*Pa)
    function dendp(mat::Material, tc::Float64, p::Float64) end
    #焓对密度偏导(J*m3/kg2)
    function dendr(mat::Material, tc::Float64, p::Float64) end
    #密度对压力偏导(s2/m2)
    function drdp(mat::Material, tc::Float64, p::Float64) end
################################################################################
end
