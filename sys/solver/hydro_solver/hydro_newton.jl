###############################################################################
# HydroNewton
###############################################################################
    mutable struct HydroNewton <: BasicSolver
        ## data
        mesh::SYSEntity.Hydro1DMesh

        ## setting
        dt::Float64
        relaxation::Float64
        max_iteration::Int64
        pre_cal::Int
        #params::Dict{Symbol, Float64}

        ## Matrix
        A::SparseMatrixCSC{Float64,Int64}
        B_n::Vector{Float64}
        B::Vector{Float64}

        HydroNewton(p::SYSEntity.Plant, circuit::SYSEntity.Circuit) =
        init(new(p.circuits[circuit], 0.0, 0.05, 100, 0, sparse(Int64[],Int64[],Float64[]), Float64[], Float64[]))
    end

    function init(s::HydroNewton)
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
        s.B_n = zeros(length(s.mesh.cvs))
        s.B = zeros(length(s.mesh.cvs))
        return s
    end

    function set_dt(s::HydroNewton, dt::Float64)
        s.dt = dt
    end

    ################################################################################
    function solve(s::HydroNewton)
        print("     solve :          $(s.mesh.circuit.name)")
        outer_construct_B(s)
        for iteration in 1:s.max_iteration
            inner_construct_AB(s)
            #display(Array(s.A))
            #display(s.B)
            #println()
            dp::Vector{Float64} = s.A \ s.B
            inner_update(s, dp)
            if convergence(s) break end
            if iteration == s.max_iteration
                println("  ERROR, $iteration  max_iteration")
                exit()
            end
        end
        outer_update(s)
    end

    function convergence(s::HydroNewton)
        return maximum(s.B) <= 1e-4
    end

    function inner_construct_AB(s::HydroNewton)
        s.B = copy(s.B_n)
        #@threads
        for cv in s.mesh.cvs
            if typeof(cv) == SYSEntity.TDV
                s.A[cv.id, cv.id] = 1.0
            else
                eii::Float64 = 0.0
                aii::Float64 = 0.0
                rh = MaterialLib.dendr(s.mesh.mat, cv[:t], cv[:p])*cv[:rho]+ cv[:enthalpy]
                for (junc, another_cv, direction) in cv.links
                    if !(typeof(junc) <: SYSEntity.TDJ)
                        imc::Float64 =  2.0 * s.dt * s.dt * junc[:area] / cv.ele[:AREA] / cv.ele[:DELTA_X] / (cv.ele[:DELTA_X] + another_cv.ele[:DELTA_X])
                        aij::Float64 = imc * (rh - junc[:enthalpy])
                        aii -= aij
                        eii += imc * (cv[:p] - another_cv[:p]) * junc[:enthalpy]
                        s.A[cv.id, another_cv.id] = aij
                    end
                end
                aii += (MaterialLib.dendp(s.mesh.mat, cv[:t], cv[:p])*cv[:rho] - 1.0)
                eii += cv[:rho]*cv[:enthalpy]-cv[:p]
                s.A[cv.id, cv.id] = aii
                s.B[cv.id] -= eii
            end
        end
    end

    function outer_construct_B(s::HydroNewton)
        for cv in s.mesh.cvs
            s.B_n[cv.id] = 0.0
            if typeof(cv) == SYSEntity.CV
                cn::Float64 = 0.0
                p::Float64 = 0.0
                for (junc, another_cv, direction) in cv.links
                    if !(typeof(junc) <: SYSEntity.TDJ)
                        p -= direction*junc[:p]*junc[:area]
                    end
                    cn -= direction*junc[:enthalpy]*junc[:rho]*junc[:area]*junc[:cv_n]
                end
                cn *= s.dt / cv.ele[:DELTA_X] / cv.ele[:AREA]
                cn += cv[:p]-cv[:rho]*cv[:enthalpy]
                cn -= p * s.dt * cv[:v] / cv.ele[:AREA]/ cv.ele[:DELTA_X]
                cn -= s.dt*cv[:q]
                #println(cv[:q],"  ",cv[:q_bak])
                #cn -= s.dt*(cv[:v]*cv[:flow_resistance])
                s.B_n[cv.id] = -cn
            end
        end
    end

    function inner_update(s::HydroNewton, delta_p::Vector{Float64})
        foreach(ii->s.mesh.cvs[ii][:p]+=delta_p[ii], 1:length(delta_p))                  # inner_update cvs' pressure
        if check_relaxation(s)
            foreach(ii->s.mesh.cvs[ii][:p]-=(1-s.relaxation)*delta_p[ii], 1:length(delta_p))
        end
        foreach(junc->inner_update_junction_velocity(junc, s.dt), s.mesh.juncs)          # inner_update junctions' velocity
        foreach(cv->inner_update_cv_property(cv, s.mesh.mat, s.dt), s.mesh.cvs)     # inner_update cvs' property
    end

    function check_relaxation(s::HydroNewton)
        if s.pre_cal > 0
            return true
        else
            for ii in 1:length(s.mesh.cvs)
                if s.mesh.cvs[ii][:p] < 0.0
                    #return true
                    return true
                end
            end
            return false
        end
    end

    function outer_update(s::HydroNewton)
        s.pre_cal -= 1
        foreach(outer_update_junction_property, s.mesh.juncs)
        foreach(cv->outer_update_cv_velocity(cv, s.mesh.mat), s.mesh.cvs)
        #foreach(outer_update, c.elements)
        foreach(junc->outer_bak_cv_n(junc, s.dt), s.mesh.juncs)
    end

    function inner_update_junction_velocity(junc::SYSEntity.Junction, dt::Float64)
        junc[:v] = junc[:cv_n] + dt*(junc.prev[:p]-junc.next[:p])/junc[:rho]/junc.ele[:DELTA_X]
        #println("junc v:", junc.id, "  ",junc[:v])
    end

    function inner_update_cv_property(cv::SYSEntity.CV, mat::MaterialLib.Material, dt::Float64)
        mv_sum::Float64 = sum(link->link[1][:rho]*link[1][:v]*link[1][:area]*link[3], cv.links)
        delta_rho::Float64 = dt * mv_sum / cv.ele[:DELTA_X] / cv.ele[:AREA]
        delta_enthalpy::Float64 = (cv[:p]-cv[:p_n])*MaterialLib.dendp(mat, cv[:t], cv[:p])+delta_rho*MaterialLib.dendr(mat, cv[:t], cv[:p])
        cv[:rho] = cv[:rho_n] + delta_rho
        cv[:enthalpy] = cv[:enthalpy_n]+ delta_enthalpy
        cv[:t] = MaterialLib.temperature(mat, cv[:enthalpy], cv[:p])
    end

    function outer_update_junction_property(junc::SYSEntity.Junction)
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

    function outer_update_cv_velocity(cv::SYSEntity.CV, mat::MaterialLib.Material)
        # bak property
        cv[:p_n] = cv[:p]
        cv[:rho_n] = cv[:rho]
        cv[:enthalpy_n] = cv[:enthalpy]
        # remove from solver later
        cv[:q] = cv[:q_bak]
        cv[:lamda] = MaterialLib.lamda(mat, cv[:t], cv[:p])
        cv[:cp] = MaterialLib.c_p(mat, cv[:t], cv[:p])
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

    function outer_update_cv_velocity(cv::SYSEntity.TDV, mat::MaterialLib.Material)
        cv[:v] = cv.links[1][1][:v]
    end

    # Update cv_bak
    function outer_bak_cv_n(junc::SYSEntity.Junction, dt::Float64)
        dz2::Float64 = junc.prev.ele[:DELTA_X]+junc.next.ele[:DELTA_X]
        tmp_a::Float64 = ((junc.prev[:v]+junc.next[:v])*(junc.prev[:v]-junc.next[:v]))/dz2
        gravity_a::Float64 = GRAVITY * sind(junc.ele[:THETA])
        flow_resistence_a::Float64 = (junc.prev[:flow_resistance]*junc.prev.ele[:DELTA_X] + junc.next[:flow_resistance]*junc.next.ele[:DELTA_X])/dz2/junc[:rho]
        junc[:cv_n]= junc[:v] + dt * (tmp_a - gravity_a - flow_resistence_a - junc[:head_source_a])
    end

    function inner_update_cv_property(cv::SYSEntity.TDV, mat::MaterialLib.Material, dt::Float64) nothing end
    function inner_update_junction_velocity(junc::SYSEntity.TDJ, dt::Float64) nothing end
    function outer_update_junction_property(junc::SYSEntity.TDJ) nothing end
    function outer_bak_cv_n(junc::SYSEntity.TDJ, dt::Float64) nothing end
