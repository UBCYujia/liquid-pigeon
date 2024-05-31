using DataFrames
using CSV

function generate_clones(num_copies::Int)
    values1 = (10, 1, 9, 2, 8, 3, 7, 4)
    values2 = (1, 10, 2, 9, 3, 8, 4, 7)

    df = DataFrame(bin_id=String[], Clone_0=Int[], Clone_1=Int[], Clone_2=Int[], Clone_3=Int[], Clone_4=Int[], Clone_5=Int[], Clone_6=Int[], Clone_7=Int[])

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

data = generate_clones(5000)

CSV.write("data/8-clones-10000.tsv", data, delim='\t')

println(first(data, 10))  
println(last(data, 10))   
