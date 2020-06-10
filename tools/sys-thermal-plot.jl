    using CSV
    using Plots

    if length(ARGS)<1
        println("Please input Setting-File names")
        exit()
    else
        input_file = string(ARGS[1])
    end

    plot_df = CSV.read(input_file)
    final = [x for x in values(last(plot_df))[2:end-1]]
    array = reshape(final, (Int(length(final)/5), 5))
    array_reverse = copy(array)
    for i in 1:size(array,1)
        array_reverse[i,:] = array[40+1-i,:]
    end
    heatmap(1:5, 1:Int(length(final)/5), array_reverse, aspect_ratio=0.3)
