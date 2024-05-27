using DataFrames
using CSV

function generate_clones(num_copies::Int)
    values1 = (10, 1, 5)
    values2 = (1, 10, 5)

    df = DataFrame(bin_id=String[], Clone_0=Int[], Clone_1=Int[], Clone_2=Int[])

    for i in 1:num_copies
        bin_label = "bin_$(i-1)"
        push!(df, (bin_label, values1...))
    end

    for i in (num_copies + 1):(2 * num_copies)
        bin_label = "bin_$(i-1)"
        push!(df, (bin_label, values2...))
    end

    return df
end

data = generate_clones(25)

CSV.write("data/3-clones-50.tsv", data, delim='\t')

println(first(data, 10))  
println(last(data, 10))   
