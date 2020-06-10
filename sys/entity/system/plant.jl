    ##------------------------------------------------------------------------------
    ## Plant
    ##------------------------------------------------------------------------------
    mutable struct Plant
        name::Symbol
        circuits::Dict{Circuit, Hydro1DMesh}
        structures::Dict{Structure, Thermal2DMesh}

        Plant() = new(:NONE, Dict{Circuit, Hydro1DMesh}(),
                      Dict{Structure, Thermal2DMesh}())
        Plant(plant_configuration::NamedTuple) = init(Plant(), plant_configuration)
    end

    function init(p::Plant, plant_configuration::NamedTuple)
        for key in keys(plant_configuration[:CIRCUIT])
            circuit = Circuit(key, plant_configuration[:CIRCUIT][key])
            push!(p.circuits, (circuit => Hydro1DMesh(circuit)))
        end
        for key in keys(plant_configuration[:THERMAL_STRUCTURE])
            structure = Structure(key, plant_configuration[:THERMAL_STRUCTURE][key])
            push!(p.structures, (structure => Thermal2DMesh(structure)))
        end
        set_structure_boundary(p)
        return p
    end

    function set_structure_boundary(p::Plant)
        for (structure, mesh) in p.structures
            # inner boundary
            if structure[:INNER_BOUNDARY_TYPE] == 4
                grids = mesh[:INNER]
                cvs = get_cvs(p, structure[:INNER_BOUNDARY_LINK])
                if cvs[1].ele[:THETA] < 0.0
                    foreach(idx->link(grids[idx], cvs[idx]), 1:length(grids))
                else
                    foreach(idx->link(grids[1+length(grids)-idx],cvs[idx]), 1:length(grids))
                end
            end
            # outer boundary
            if structure[:OUTER_BOUNDARY_TYPE] == 4
                grids = mesh[:OUTER]
                cvs = get_cvs(p, structure[:OUTER_BOUNDARY_LINK])
                if cvs[1].ele[:THETA] < 0.0
                    foreach(idx->link(grids[idx], cvs[idx]), 1:length(grids))
                else
                    foreach(idx->link(grids[1+length(grids)-idx],cvs[idx]), 1:length(grids))
                end
            end
        end
    end

    function get_element(p::Plant, ele::Symbol)
        for circuit in keys(p.circuits)
            if ele in keys(circuit)
                return circuit[ele]
            end
        end
        nothing
    end

    function get_cvs(p::Plant, ele::Symbol)
        for circuit in keys(p.circuits)
            if ele in keys(circuit)
                return p.circuits[circuit][ele]
            end
        end
        nothing
    end

    function get_all_cvs(p::Plant)
        cvs = Vector{SYSEntity.Node}()
        for circuit in keys(p.circuits)
            cvs = [cvs; p.circuits[circuit][]]
        end
        return cvs
    end

    function get_structure(p::Plant, key::Symbol)
        for structure in keys(p.structures)
            if key == p.name
                return structure
            end
        end
        nothing
    end

    function get_grids(p::Plant, key::Symbol)
        for structure in keys(p.structures)
            if key == p.name
                return p.structures[structure][]
            end
        end
        nothing
    end

    function get_all_grids(p::Plant)
        grids = Vector{SYSEntity.Grid}()
        for structure in keys(p.structures)
            grids = [grids; p.structures[structure][]]
        end
        return grids
    end

    function display_info(p::Plant)
        for (cir,mesh) in p.circuits
            println("    $(cir.name),  Element Num: $(length(cir.elements))")
        end
        for (structure,grids) in p.structures
            println("    $(structure.name)")
        end
        println("***********************************************************************************")
    end
