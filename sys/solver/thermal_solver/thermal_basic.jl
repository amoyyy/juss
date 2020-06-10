###############################################################################
# ThermalBasic
###############################################################################
    mutable struct ThermalBasic <: BasicSolver
        ## data
        mesh::SYSEntity.Thermal2DMesh

        ## settings
        dt::Float64
        max_iteration::Int64

        ## matrix
        A::SparseMatrixCSC{Float64,Int64}
        B::Vector{Float64}
        ThermalBasic(p::SYSEntity.Plant, structure::SYSEntity.Structure) =
        init(new(p.structures[structure],
                 0.0, 30, sparse(Int64[],Int64[],Float64[]), Float64[]))
    end

    function init(s::ThermalBasic)
        ## Matrix
        s.A = sparse(1:length(s.mesh.grids), 1:length(s.mesh.grids), rand(length(s.mesh.grids)))
        s.B = zeros(length(s.mesh.grids))
        #
        return s
    end

    function set_dt(s::ThermalBasic, dt::Float64)
        s.dt = dt
    end

    function solve(s::ThermalBasic)
        print("     solve :          $(s.mesh.structure.name)")
        fetch_hydro(s)
        for iteration in 1:s.max_iteration
            construct_AB(s)
            T::Vector{Float64} = s.A \ s.B
            update(s, T)
            if convergence(s) break end
        end
        feedback_hydro(s)
    end

    function convergence(s::ThermalBasic)
        return maximum(s.B) <= 1e-4
    end

    function construct_AB(s::ThermalBasic)
        for g1 in s.mesh.grids
            a_p::Float64 = 0.0
            b::Float64 = 0.0
            for link in g1.links
                g2 = link[1]
                f1 = link[2]
                r::Float64 = g1[:cordin2]+f1*link[3]*0.5*g1[:delta_cordin2]
                dr::Float64 = f1*g1[:delta_cordin1]+(1-f1)*g1[:delta_cordin2]
                hr::Float64 = f1*(g1[:heat_resistence2]+g2[:heat_resistence2])+(1-f1)*(g1[:heat_resistence1]+g2[:heat_resistence1])
                coeff::Float64 = 2.0*r*dr/hr
                a_p += coeff
                s.A[g1.id, g2.id] = -coeff
            end
            tmp::Float64 = g1[:rho]*g1[:cp]/s.dt
            a_p += (tmp-g1[:q_sp])*g1[:delta_v]
            b = (tmp*g1[:t]+g1[:q_sc])*g1[:delta_v]
            # Boundary handling
            if g1.type != 0
                A::Float64 = g1[:delta_cordin1]*g1[:boundary_cordin]
                if g1.type == 1
                    b += 2.0*(g1[:twall]-g1[:t])*A/g1[:heat_resistence2]
                elseif g1.type == 2
                    b += g1[:qwall]*A
                elseif g1.type == 3 || g1.type == 4
                    h_modify::Float64 = 1.0/(1.0/g1[:h_coeff]+0.5*g1[:heat_resistence2])
                    b += g1[:tfluid]*h_modify*A
                    a_p += h_modify*A
                end
            end
            # coeff output
            s.A[g1.id, g1.id] = a_p
            s.B[g1.id] = b
        end
    end

    function update(s::ThermalBasic, T::Vector{Float64})
        foreach(ii->s.mesh.grids[ii][:t]=T[ii], 1:length(T))
        foreach(grid->update_property(grid, s.mesh.mat), s.mesh.grids)
    end

    function fetch_hydro(s::ThermalBasic)
        foreach(SYSEntity.fetch_boundary, s.mesh.grids)
    end

    function feedback_hydro(s::ThermalBasic)
        foreach(SYSEntity.feedback_boundary, s.mesh.grids)
    end

    function update_property(g::SYSEntity.Grid, mat::MaterialLib.Material)
        g[:rho] = MaterialLib.ro(mat, g[:t], 0.0)
        g[:lamda] = MaterialLib.lamda(mat, g[:t], 0.0)
        g[:cp] = MaterialLib.c_p(mat, g[:t], 0.0)
        g[:heat_resistence1] = g[:delta_cordin1] / g[:lamda]
        g[:heat_resistence2] = g[:delta_cordin2] / g[:lamda]
    end
