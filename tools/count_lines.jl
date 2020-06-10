root_dir = "./"
total_files = 0
total_lines = 0
for (root, dirs, files) in walkdir(root_dir)
    for file in files
        if (occursin(".jl", file))
            lines::Int = 0
            global total_files += 1
            open(joinpath(root, file), "r") do io
                lines = countlines(io)
                global total_lines += lines
            end
            println(total_files," ",joinpath(root, file),"  ", lines)
        end
    end
end

println("Num_of_files: ", total_files,"    Num_of_lines: ", total_lines)
