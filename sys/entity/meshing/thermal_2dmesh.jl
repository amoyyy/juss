    mutable struct Thermal2DMesh <: Mesh
        structure::Structure
        mat::MaterialLib.Material
        grids::Matrix{Node}
        #
        Thermal2DMesh(s::Structure) =
        init(new(s, MaterialLib.Steel(), Matrix{Node}(undef, size(s))))
    end

    function init(m::Thermal2DMesh)
        m.mat = MaterialLib.getMat(m.structure[:MATERIAL_MEDIA])
        (axial_num, radial_num) = size(m.grids)
        # add nodes
        id_idx::Int64 = 1
        for ii in 1:radial_num
            for jj in 1:axial_num
                set_node(m, Grid(id_idx, m.structure, jj, ii))
                id_idx += 1
            end
        end
        # link nodes
        for ii in 1:radial_num
            for jj in 1:axial_num-1
                link(m[jj, ii], m[jj+1, ii])
            end
        end
        for jj in 1:axial_num
            for ii in 1:radial_num-1
                link(m[jj, ii], m[jj, ii+1])
            end
        end
        foreach(grid->init(grid, m.mat), m.grids)
        return m
    end

    function set_node(m::Thermal2DMesh, node::Grid)
        m.grids[node.id_axial, node.id_radial] = node
    end

    function Base.getindex(m::Thermal2DMesh, axial::Int64, radial::Int64)
        return m.grids[axial, radial]
    end

    function Base.getindex(m::Thermal2DMesh, cordin::Symbol)
        if cordin == :INNER
            return m.grids[:,1]
        elseif cordin == :OUTER
            return m.grids[:,end]
        end
    end

    function Base.getindex(m::Thermal2DMesh, id::Int64)
        return m.grids[radial]
    end

    function Base.getindex(m::Thermal2DMesh)
        return m.grids[1:end]
    end
