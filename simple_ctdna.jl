using BridgeStan
using Pigeons
using Random
using CSV
using DataFrames
using JSON
using InferenceReport
struct StanCtDNAModel end

# initialize Stan model with ctdna data
function init_stan_ctdna(n, num_clones, clone_cn_profiles, ctdna, scale)
    stan_file = "model/simple_model.stan"
    stan_data = Pigeons.json(; n=n, num_clones=num_clones, clone_cn_profiles=clone_cn_profiles, ctdna=ctdna, scale=scale)

    return StanLogPotential(stan_file, stan_data, StanCtDNAModel())
end

ctdna_path = "data/ctdna.tsv"
clones_path = "data/2-clones-simple.tsv"
stan_file = "model/simple_model.stan"

ctdna_data = CSV.read(ctdna_path, DataFrame, delim='\t',header=false)
clones_data = CSV.read(clones_path, DataFrame, delim='\t')
# check the first few rows
println(first(ctdna_data, 5))
println(first(clones_data,5))



n = size(clones_data)[1]
num_clones = size(clones_data)[2]-1
scale = 1.0 #for simplicity for now
# trigger error if using Pigeons.json()
data_dict = Dict(
    "n" => n,
    "num_clones" => num_clones,
    "clone_cn_profiles" => transpose(Matrix(clones_data[:, 2:end])),  # Transpose to fit matrix[n, num_clones]
    "ctdna" => Vector(ctdna_data[:,1]),  # This assumes you've added a 'ctdna' column to your DataFrame
    "scale" => 1.0
)
json_data = JSON.json(data_dict)

pt = pigeons(target = StanLogPotential(stan_file,json_data))

report(pt)