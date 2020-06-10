SRC_ROOT_DIR = "/home/amos/Projects/Juss_2.0"

function ADD_LOAD_PATH(ROOT_DIR::String)
    MODULE_PATHS = String[]
    if isdir(ROOT_DIR)
        push!(MODULE_PATHS, ROOT_DIR)
        for (root, dirs, files) in walkdir(ROOT_DIR)
            for dir in dirs
                ADD_DIR = joinpath(root,dir)
                if !occursin("cases", ADD_DIR)
                    push!(MODULE_PATHS, ADD_DIR)
                end
            end
        end
    end
    foreach(path->if !(path in LOAD_PATH)
                      push!(LOAD_PATH,path)
                  end,
                  MODULE_PATHS)
end

ADD_LOAD_PATH(SRC_ROOT_DIR)
