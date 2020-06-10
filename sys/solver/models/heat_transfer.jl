###############################################################################
## HeatTransfer
###############################################################################
export heat_transfer_coeff
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
    function heat_transfer_coeff(rho::Float64, lamda::Float64, cp::Float64, v::Float64, d::Float64)
        #  Dcs=4*Acs/Pst;
        #  Pe=abs(Wcs)*Ccs*Dcs/(Kcs*Acs)
        #  Nu=5.0+0.025*Pe**0.8
        #  Hst=Kcs/Dcs*Nu !流体与壁面传热系数
        #  Hst=1/(1/Hst+Dtu/(2*Ktt)) !流体与管壁的传热系数（修正后）
        alpha::Float64 = lamda/rho/cp
        pe::Float64 = Pe(abs(v), d, alpha)
        nu::Float64 = Nu(pe,0.0,0.0)
        return H(nu, d, lamda)
    end

    #Nu = (pe::Float64, c1::Float64,  c2::Float64) -> 5.0+0.025*(pe^0.8)
    #Nu = (pe::Float64, c1::Float64,  c2::Float64) -> 5.3+0.01*(pe^0.85)
    function Nu(pe::Float64, c1::Float64,  c2::Float64)
        return 5.0+0.025*(pe^0.8)
    end

    #Pe = (v::Float64, l::Float64,  a::Float64) ->  v*l/a
    function Pe(v::Float64, l::Float64,  a::Float64)
        return v*l/a
    end

    #H = (nu::Float64, x::Float64,  k::Float64) ->  k*nu/x
    function H(nu::Float64, x::Float64,  k::Float64)
        return k*nu/x
    end

    """//H20对流换热系数
    //cv 是流体信息， v[0]: wall temperature, v[1] heat flux between cv and tv. just for water
    double water::heat_exchange_coeff(control_volume *cv, std::vector<double> v) {
        //全局变量（需要从外界获得的变量）
        double p = cv->_pressure;         //压力
        double t = cv->_temperature;             //温度
        double dens = cv->_density;      //密度
        double Diameter =
                4 * cv->_area / cv->_ele->_structure_param[structure_param_index::WETTED_PERIMETER_INDEX];     //水力直径
        double Velocity = std::fabs(cv->_velocity);     //速度
        double myq = std::fabs(v[1]);     //每个节点的热流密度
        double xfraction = vap_fraction(cv->_enthalpy, cv->_pressure);     //含气率 范围是0~1
        double my_twall = v[0];    //壁面温度
    //局部变量（程序内部使用的中间变量）
        double temp_Visc = VISC(t, dens);
        double Hcon = THCOND(t, dens);
        double Re_number = (Diameter / temp_Visc) * dens * Velocity;
        double temp_G = dens * Velocity;
        double chfq = chf(xfraction, temp_G, cv->_pressure);  //计算得到的临界热流密度
        double htransfer = 0.0;
    //进行计算
        if (xfraction <= 1.e-30) {
            //过冷段   D-B公式
            double Pr_number = PRC(p, t);
            double Nu_number = 0.023 * pow(Re_number, 0.8) * pow(Pr_number, 0.4);
            htransfer = Nu_number / Diameter * Hcon;
            //std::cout << "1: " << cv->_index_cv_of_element <<" : "<< htransfer << std::endl;
        } else if (xfraction > 1.e-30 && xfraction < 1.0) {
            if (myq < chfq) {
                //表面热流密度低于临界热流密度，处于泡核沸腾，使用chen
                double tempx = pow((xfraction / (1 - xfraction)), 0.9) * pow((SGSV(p) / SFSV(p)), 0.5) *
                               pow((VISSG(p) / VISCL(p, 1)), 0.1);
                double tempF;
                if (tempx <= 0.1) {
                    tempF = 1.0;
                } else {
                    tempF = 2.35 * pow((tempx + 0.213), 0.736);
                }
                double tempRetp = temp_G * (1 - xfraction) * Diameter / VISCL(p, 1) * pow(tempF, 1.25) / 10000;
                double tempS;
                if (tempRetp < 32.5) {
                    tempS = 1 / (1 + 0.12 * pow(tempRetp, 1.14));
                } else if (tempRetp >= 32.5 && tempRetp < 70.0) {
                    tempS = 1 / (1 + 0.42 * pow(tempRetp, 0.78));
                } else {
                    tempS = 0.1;
                }
                double tempCp = (HLPT(p, t - 2.) - HLPT(p, t - 4.0)) / 2.0;
                double temp_rhof = 1 / SFSV(p);
                double temp_hmac = 0.023 * tempF * pow(THCOND(t, temp_rhof), 0.6) * pow(temp_G, 0.8) *
                                   pow((1. - xfraction), 0.8) * pow(tempCp, 0.4) /
                                   (pow(VISCL(p, 1), 0.4) * pow(Diameter, 0.2));
                double temp_hmic = 0.00122 * tempS * pow((THCOND(t, temp_rhof)), 0.79) * pow(tempCp, 0.45) *
                                   pow(temp_rhof, 0.49) / (pow(surftense(t), 0.5) * (pow((VISCL(p, 1)), 0.29) *
                                                                                     pow(((SGHP(p) - SFHP(p))), 0.24) *
                                                                                     pow((1 / SGSV(p)), 0.24))) *
                                   pow((abs(my_twall - TSATP(p))), 0.24) *
                                   //pow((abs(PSATT(min(my_twall, t + 15.0)) - p)), 0.75);
                                   pow((abs(PSATT(my_twall) - p)), 0.75);
                htransfer = temp_hmac + temp_hmic;
                //ETEC correlations
    //            htransfer = pow(myq, 0.5) * exp(p / 8.69e6) / 0.02253;
    //            std::cout << "2: " << cv->_index_cv_of_element <<" : "<< htransfer << std::endl;
            } else {
                //表面热流密度大于临界热流密度，处于膜态沸腾，使用Groneveld公式
                double temp_YY = 1.0 - 0.1 * pow((1 - xfraction), 0.4) * pow((SGSV(p) / SFSV(p) - 1), 0.4);
                temp_YY = max(temp_YY, 0.1);
                double cp_g = (HGPT(p, TSATP(p) + 5.) - HGPT(p, TSATP(p) + 3.)) / 2.;
                //Groneveld公式
                double Nu_number = 0.052 * pow((temp_G * Diameter / (VISSG(p)) *
                                                (xfraction + (1 - xfraction) * SFSV(p) / SGSV(p))), 0.688) *
                                   pow((cp_g * (VISSG(p)) / Hcon), 1.26) * pow(temp_YY, -1.06);
                //ETEC correlations
    //            Nu_number = 0.0193 * pow(temp_G * Diameter / VISSG(p), 0.8) * pow(PRG(p), 1.23) *
    //                        pow(xfraction + (1 - xfraction) * SFSV(p) / SGSV(p), 0.68) * pow(SFSV(p) / SGSV(p), 0.068);
                htransfer = Nu_number / Diameter * Hcon;
    //            std::cout << "3: " << cv->_index_cv_of_element <<" : "<< htransfer << std::endl;
            }
        } else {
            //过热汽区域   Sider-Tate公式
            double Pr_number = PRC(p, t);
            double Nu_number = 0.023 * pow(Re_number, 0.8) * pow(Pr_number, 0.33) *
                               pow((temp_Visc / (VISC(my_twall, dens))), 0.14);
            //ETEC correlations
            //Nu_number = 0.0073 * pow(Re_number, 0.886) * pow(Pr_number, 0.61);
            //std::cout << "Pr:" << Pr_number << " " << v[0] << std::endl;
            htransfer = Nu_number / Diameter * Hcon;
            //std::cout << "4: " << cv->_index_cv_of_element <<" : "<< htransfer << std::endl;
        }
        return htransfer;
    }
"""
