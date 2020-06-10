    struct Sodium <: Material end
    ##------------------------------------------------------------------------------
    ## Sodium
    ##------------------------------------------------------------------------------
    function enth(mat::Sodium, t::Float64, p::Float64)#焓(J/kg)
        tf=32.0+1.8*t
        tr=tf+459.67
        patm=p/101325.0
        if tr>2059.67
            psat=6881760.2/(tr^0.61344)*exp(-22981.96/tr)
        else
            psat=3032660.0/sqrt(tr)*exp(-23073.3/tr)
        end
        hls=-29.02+tr*(0.389352+tr*(-0.552955e-4+tr*0.113726e-7))
        d=ro(mat,t,p)/16.01846
        dlp=(1.8105e-10*tf-5.744e-7)*tf-7.9504e-3
        hlc=hls+1.0/d*(1.0+tr/d*dlp)*(patm-psat)*2.721308
        return hlc*2326.0
    end

    function ro(mat::Sodium, tc::Float64, p::Float64)#密度(kg/m3)
        return 9.50076E2+tc*(-.22976+tc*(-1.46049e-5+5.63788E-9*tc))
    end

    function c_p(mat::Sodium, tc::Float64, p::Float64)#比热(J/kg*K)
        tf=32.0+1.8*tc
        tr=tf+459.67
        patm=p/101325.0
        if tr>2059.67
            psat=6881760.2/(tr^0.61344)*exp(-22981.96/tr)
            dpsadt=psat/tr*1.8*(-0.61344+22981.96/tr)
        else
            psat=3032660.0/sqrt(tr)*exp(-23073.3/tr)
            dpsadt=psat/tr*1.8*(-0.5+23073.3/tr)
        end
        hls=-29.02+tr*(0.389352+tr*(-0.552955e-4+tr*0.113726e-7))
        dhlsdt=((3.411178e-8*tr-1.105991e-4)*tr+0.389352)*1.8
        rol=9.50076e2+tc*(-2.2976e-1+tc*(-1.46049e-5+5.63788e-9*tc))
        d=rol/16.01846
        dsq=d*d
        ddt=(-2.2976e-1+tc*(-2.92098e-5+1.69136e-8*tc))/16.01846
        dlp=(1.8105e-10*tf-5.744e-7)*tf-7.9504e-3
        ddlpdt=(3.621e-10*tf-5.7744e-7)*1.8
        f1=1.0+tr/d*dlp
        f2=(patm-psat)/d*2.721308
        df1dt=-tr*dlp*ddt/dsq+tr/d*ddlpdt+dlp/d*1.8
        df2dt=-(dpsadt/d+(patm-psat)*ddt/dsq)*2.721308
        dhlcdt=dhlsdt+df1dt*f2+df2dt*f1
        return dhlcdt*2326.0
    end

    function temperature(mat::Sodium, enthalpy::Float64, pressure::Float64)
        #result = -290.3879954627355 + 0.0007806716530000233*enthalpy - 7.304009366423947e-07*pressure
        #result = -291.6627868932519 + 0.0007823376026054913*enthalpy - 7.035246803038797e-07*pressure
        #result = -291.6627868932519 + 0.0007823376026054913*enthalpy
        # 2018-10-11 MOX 燃料组件调试
        result = 0.0007860048*enthalpy-298.2
        return result
    end

    function dr(mat::Sodium, temp::Float64)
        return -0.22976-2.92098e-5*temp+1.691364e-8*temp*temp
    end

    function tpsat(mat::Sodium, pressure::Float64)#求钠的饱和温度 temp(degC)
        pa=pressure*0.986924/1e5
        #println(pa,"  ",pressure)
        a=log10(pa)
        a1=a-6.4818
        a2=a-6.8377
        a3=a-1.36041
        ttest=2e3
        tsatp=0
        for i in 1:59
            tsatp=abs(ttest)
            if ttest<=2059.7
                ttest=-10020.6/(a1+0.5*log10(ttest))
            elseif ttest<=2959.7
                ttest=-9980.94/(a2+0.61344*log10(ttest))
            else
                ttest=-8178.27/(a3-0.789*log10(ttest))
            end
            if abs((ttest-tsatp)/ttest)<=1e-4
                break
            end
        end
        return (tsatp-459.7-32)/1.8
    end

    function eta(mat::Sodium, temp::Float64, pressure::Float64)#动力粘度(kg/m*s)
        tr=temp*1.8+32+459.7
        tf=temp*1.8+32
        if temp<=97.8#固体
            result = -1
        elseif temp<tpsat(mat,pressure)
            result = 4.1338e-4*exp(2.302585*(1.0203+397.17/tr-0.4925*log10(tr)))
        else
            result = 4.134e-4*(0.03427+8.176e-6*tf)
        end
        return result
    end

    function lamda(mat::Sodium, temp::Float64, pressure::Float64)#导热系数(W/m*K)
        tf=temp*1.8+32
        result = 0.0
        if temp<=97.8   #固体
            result=(1.356-0.00167*temp)*100
        elseif temp<tpsat(mat,pressure)
        #elseif temp<882+273.15
            result=(54.306-1.878e-2*tf+2.0914e-6*tf*tf)*1.7307
        else
            result=(0.1639e-2+0.3977e-4*tf-0.9697e-8*tf*tf)*1.7307
        end
        return result
    end

    function dendp(mat::Sodium, tc::Float64, p::Float64)
        rol=9.50076e2+tc*(-2.2976e-1+tc*(-1.46049e-5+5.63788e-9*tc))
        drdt=-2.2976e-1+tc*(-2.92098e-5+1.69136e-8*tc)
        alfa=-drdt/rol
        return (1.0-(tc+273.15)*alfa)/rol
    end

    function dendr(mat::Sodium, tc::Float64, p::Float64)
        soldrh=c_p(mat,tc,p)/(-2.2976e-1+tc*(-2.92099e-5+1.69136e-8*tc))
        return soldrh
    end

    function drdp(mat::Sodium, tc::Float64, p::Float64)
        rol=9.50076e2+tc*(-2.2976e-1+tc*(-1.46049e-5+5.63788e-9*tc))
        cpl=c_p(mat,tc,p)
        drdt=-2.2976e-1+tc*(-2.92098e-5+1.69136e-8*tc)
        alfa=-drdt/rol
        return alfa/cpl*(1.0-alfa*(tc+273.15))
    end

    function cs(mat::Sodium, tc::Float64, p::Float64)
        return 1.0/sqrt(drdp(mat, tc, p))
    end
