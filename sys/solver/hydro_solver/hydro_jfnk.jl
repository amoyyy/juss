###############################################################################
# HydroJFNK
###############################################################################
    mutable struct HydroJFNK <: BasicSolver
        ## data
        mesh::SYSEntity.Hydro1DMesh

        # #settings
        dt::Float64
        relaxation::Float64
        max_iteration::Int64
        pre_time::Int
        #params::Dict{Symbol, Float64}

        ## matrix
        A::SparseMatrixCSC{Float64,Int64}
        B::Vector{Float64}
        HydroJFNK(p::SYSEntity.Plant, circuit::SYSEntity.Circuit) =
        init(new(p.circuits[circuit], 0.0, 1.0, 9, 5, sparse(Int64[],Int64[],Float64[]), Float64[]))
    end

    function init(s::HydroJFNK)
        ii::Vector{Int64} = Int64[]
        jj::Vector{Int64} = Int64[]
        for cv in s.mesh.cvs
            if typeof(cv) == SYSEntity.TDV
                push!(ii, cv.id)
                push!(jj, cv.id)
            else
                for (junc, another_cv, direction) in cv.links
                    if !(typeof(junc) <: SYSEntity.TDJ)
                        push!(ii, cv.id)
                        push!(jj, another_cv.id)
                    end
                end
            end
        end
        vv::Vector{Float64} = zeros(length(ii))
        s.A = sparse(ii, jj, vv)
        s.B = zeros(size(s.A, 1))
        return s
    end

    function set_dt(s::HydroJFNK, dt::Float64)
        s.dt = dt
    end

    ################################################################################
    function solve(s::HydroJFNK)
        print("     solve :          $(s.mesh.circuit.name)")
        k_divide = get_relaxation_coeff(s)
        for i in 1:k_divide
            construct_AB(s, 1.0/k_divide)
            D  = Array(diag(s.A))
            M2 = x -> D.\x
            x = KrylovMethods.gmres(s.A, s.B, 5, tol=1e-7, maxIter=200, M=M2)
            update_p(s, x[1])
            update(s)
        end
        restore_cv_property(s)
    end

    function construct_AB(s::HydroJFNK, coeff::Float64)
        #@threads
        for cv in s.mesh.cvs
            if typeof(cv) == SYSEntity.TDV
                s.A[cv.id, cv.id] = 1.0
                s.B[cv.id] = cv[:p]
            else
                r_p::Float64 = MaterialLib.drdp(s.mesh.mat, cv[:t], cv[:p])
                en_p::Float64 = MaterialLib.dendp(s.mesh.mat, cv[:t], cv[:p])
                en_r::Float64 = MaterialLib.dendr(s.mesh.mat, cv[:t], cv[:p])
                aii::Float64 = (en_r * cv[:rho] + cv[:enthalpy]) * r_p + cv[:rho] * en_p - 1.0
                eii::Float64 = aii * cv[:p]
                tmp::Float64 = 0.0
                pa::Float64 = 0.0
                for (junc, another_cv, direction) in cv.links
                    if !(typeof(junc) <: SYSEntity.TDJ)
                        imc::Float64 =  2.0 * s.dt * s.dt * junc[:area] / cv.ele[:AREA] / cv.ele[:DELTA_X] / (cv.ele[:DELTA_X] + another_cv.ele[:DELTA_X])
                        aij::Float64 = - imc * junc[:enthalpy]
                        aii -= aij
                        pa -= direction*junc[:area] * junc[:p] * cv[:v]
                        s.A[cv.id, another_cv.id] = aij
                    end
                    tmp += direction*junc[:area]*junc[:enthalpy]*junc[:rho]*junc[:cv_n]
                end
                eii += (tmp + pa) * s.dt / cv.ele[:DELTA_X] / cv.ele[:AREA]
                eii += s.dt * cv[:q] * coeff
                #eii += s.dt*(cv[:v]*cv[:flow_resistance])
                s.A[cv.id, cv.id] = aii
                s.B[cv.id] = eii
            end
        end
    end

    function get_relaxation_coeff(s::HydroJFNK)
        max_heatflux = 1e8
        max_q = maximum([abs(cv[:q]) for cv in s.mesh.cvs])
        if s.pre_time>0
            s.pre_time -= 1
            return 1+floor(max_q/max_heatflux)
        else
            if max_q > max_heatflux*3
                println("      ",max_q)
                return 1+floor(max_q/max_heatflux)
            end
            return 1
        end
    end

    function restore_cv_property(s::HydroJFNK)
        foreach(cv->cv[:q] = cv[:q_bak], s.mesh.cvs)
    end

    function update_p(s::HydroJFNK, p::Vector{Float64}, coeff::Float64=1.0)
        for ii in 1:length(p)
            cv = s.mesh.cvs[ii]
            if typeof(cv) != SYSEntity.TDV
                cv[:p_n] = cv[:p]
                cv[:p] = coeff*p[ii]+(1-coeff)*cv[:p_n]
                #println("$(cv.id)  $(cv[:p_n])    $(cv[:p])    $(p[ii])      $(cv[:q])   $(cv[:q_bak])")
            end
        end
    end

    function update(s::HydroJFNK)
        foreach(junc->update_junction_velocity(junc, s.dt), s.mesh.juncs)
        foreach(cv->update_cv_property(cv, s.mesh.mat, s.dt), s.mesh.cvs)
        foreach(junc->update_junction_property(junc), s.mesh.juncs)
        foreach(cv->update_cv_velocity(cv), s.mesh.cvs)
        foreach(junc->bak_cv_n(junc, s.dt), s.mesh.juncs)
    end

    function update_junction_velocity(junc::SYSEntity.Junction, dt::Float64)
        junc[:v] = junc[:cv_n] + dt*(junc.prev[:p]-junc.next[:p])/junc[:rho]/junc.ele[:DELTA_X]
    end

    function update_cv_property(cv::SYSEntity.CV, mat::MaterialLib.Material, dt::Float64)
        mv_sum::Float64 = sum(link->link[1][:rho]*link[1][:v]*link[1][:area]*link[3], cv.links)
        delta_rho::Float64 = dt * mv_sum / cv.ele[:DELTA_X] / cv.ele[:AREA]
        delta_enthalpy::Float64 = (cv[:p]-cv[:p_n])*MaterialLib.dendp(mat, cv[:t], cv[:p])+delta_rho*MaterialLib.dendr(mat, cv[:t], cv[:p])
        cv[:rho] += delta_rho
        cv[:enthalpy] += delta_enthalpy
        cv[:t] = MaterialLib.temperature(mat, cv[:enthalpy], cv[:p])
        #println("$(cv[:rho])   $delta_rho    $(mv_sum),   $(cv[:enthalpy])   $delta_enthalpy   ,   $(cv[:t])")
        # remove from solver later
        cv[:lamda] = MaterialLib.lamda(mat, cv[:t], cv[:p])
        cv[:cp] = MaterialLib.c_p(mat, cv[:t], cv[:p])
    end

    function update_junction_property(junc::SYSEntity.Face)
        if junc[:v]>=0
            junc[:rho] = junc.prev[:rho]
            junc[:enthalpy] = junc.prev[:enthalpy]
            junc[:p] = junc.prev[:p]
        else
            junc[:rho] = junc.next[:rho]
            junc[:enthalpy] = junc.next[:enthalpy]
            junc[:p] = junc.next[:p]
        end
    end

    function update_cv_velocity(cv::SYSEntity.CV)
        # get weight velocity
        rhoa_sum::Float64 = 0.0
        rhoa_sum += sum(link->link[1][:rho]*link[1][:area], cv.links)

        area_sum::Float64 = 0.0
        rhoav_sum::Float64 = 0.0
        area_sum2::Float64 = 0.0
        rhoav_sum2::Float64 = 0.0
        for (junc, another_cv, direction) in cv.links
            if direction>0.0
                area_sum += junc[:area]
                rhoav_sum += junc[:rho]*junc[:area]*junc[:v]
            else
                area_sum2 += junc[:area]
                rhoav_sum2 += junc[:rho]*junc[:area]*junc[:v]
            end
        end
        cv[:v] = (rhoav_sum*area_sum+rhoav_sum2*area_sum2) / rhoa_sum / cv.ele[:AREA]
        #cv[:flow_resistance] = CalModel.s_flow_resistence_loss(cv[:v],cv[:rho],cv.ele[:HYDRO_DIAMETER],cv[:eta],0.0)
        cv[:h_coeff] = heat_transfer_coeff(cv[:rho], cv[:lamda], cv[:cp], cv[:v], cv.ele[:HYDRO_DIAMETER])
    end

    function update_cv_velocity(cv::SYSEntity.TDV)
        cv[:v] = cv.links[1][1][:v]
    end

    # Update cv_bak
    function bak_cv_n(junc::SYSEntity.Junction, dt::Float64)
        dz2::Float64 = junc.prev.ele[:DELTA_X]+junc.next.ele[:DELTA_X]
        tmp_a::Float64 = ((junc.prev[:v]+junc.next[:v])*(junc.prev[:v]-junc.next[:v]))/dz2
        gravity_a::Float64 = GRAVITY * sind(junc.ele[:THETA])
        flow_resistence_a::Float64 = (junc.prev[:flow_resistance]*junc.prev.ele[:DELTA_X] + junc.next[:flow_resistance]*junc.next.ele[:DELTA_X])/dz2/junc[:rho]
        junc[:cv_n]= junc[:v] + dt * (tmp_a - gravity_a - flow_resistence_a - junc[:head_source_a])
    end

    function update_cv_property(cv::SYSEntity.TDV, mat::MaterialLib.Material, dt::Float64) nothing end
    function update_junction_velocity(junc::SYSEntity.TDJ, dt::Float64) nothing end
    function bak_cv_n(junc::SYSEntity.TDJ, dt::Float64) nothing end
