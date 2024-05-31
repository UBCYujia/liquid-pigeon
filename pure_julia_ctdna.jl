using Pigeons
using Random
using Statistics
using Distributions
using CSV
using DataFrames
using InferenceReport

struct CtDNALogPotential
    ctdna::Vector{Float64}
    clone_cn_profiles::Matrix{Float64}
    num_clones::Int
    n::Int
    scale::Float64
end

function (log_potential::CtDNALogPotential)(params)
    start_time = time_ns()
    rho = params
    if any(x -> x < 0 || x > 1, rho) || abs(sum(rho)  - 1) > 1e-5
    # if any(x -> x < 0 || x > 1, rho) || sum(rho) != 1
        return -Inf  #ensure rho is valid
    end

    total_sum = log_potential.clone_cn_profiles * rho
    mean_total_sum = mean(total_sum)

    mu = log.(total_sum) .- log(mean_total_sum)

    log_likelihood = 0.0
    degrees_of_freedom = 2
    # dist = TDist(degrees_of_freedom)
    
    # for i in 1:log_potential.n
    #     dist_i = TDist(degrees_of_freedom) * log_potential.scale + mu[i]
    #     log_potential.ctdna[i] ~ dist_i
    #     log_likelihood += logpdf(dist_i, ctdna[i])
    # end

    dist = TDist(degrees_of_freedom)

    for i in 1:log_potential.n
        dist_i = LocationScale(mu[i], log_potential.scale, dist)
        log_likelihood += logpdf(dist_i, log_potential.ctdna[i])
    end

    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 1e9
    println("Execution of log_potential took $elapsed_time seconds")

    return log_likelihood
end


function Pigeons.initialization(log_potential::CtDNALogPotential, rng::AbstractRNG, ::Int)
    alpha = 1.0  
    rho = rand(rng, Dirichlet(log_potential.num_clones, alpha))

    @assert abs(sum(rho) - 1) < 1e-5 "density not 1!"

    return rho
end

function Pigeons.sample_iid!(log_potential::CtDNALogPotential, replica, shared)
    rng = replica.rng
    new_state = rand(rng, Dirichlet(log_potential.num_clones, 1.0))

    @assert abs(sum(new_state) - 1) < 1e-5  "density not 1!"

    replica.state = new_state
end


function load_data(ctdna_path, clones_path)
    ctdna_data = CSV.read(ctdna_path, DataFrame, delim='\t', header=false,types=[Float64])
    clones_data = CSV.read(clones_path, DataFrame, delim='\t')
    return ctdna_data, clones_data
end

function default_reference(log_potential::CtDNALogPotential)
    neutral_ctdna = randn(log_potential.n)
    neutral_cn_profiles = abs.(randn(size(log_potential.clone_cn_profiles)))
    return CtDNALogPotential(neutral_ctdna, neutral_cn_profiles, log_potential.num_clones, log_potential.n, log_potential.scale)
end

# function main()
#     ctdna_path = "data/ctdna-5000.tsv"
#     clones_path = "data/2-clones-5000.tsv"
#     ctdna_data, clones_data = load_data(ctdna_path, clones_path)

#     n = size(clones_data, 1)
#     num_clones = size(clones_data, 2) - 1
#     clone_cn_profiles = Matrix(clones_data[:, 2:end])
#     ctdna = Vector{Float64}(ctdna_data[:, 1])
#     scale = 1.0  #simplified scale

#     log_potential = CtDNALogPotential(ctdna, clone_cn_profiles, num_clones, n, scale)
#     reference_potential = default_reference(log_potential)

#     pt = pigeons(
#         target = log_potential,
#         reference = reference_potential,
#         record = [traces; record_default()]
#     )
#     # report(pt)

#     println("Model run complete.")
# end

# @time main()
function main(ctdna_paths, clones_paths)
    for (ctdna_path, clones_path) in zip(ctdna_paths, clones_paths)
        println("processing: $ctdna_path and $clones_path")
        ctdna_data, clones_data = load_data(ctdna_path, clones_path)

        n = size(clones_data, 1)
        num_clones = size(clones_data, 2) - 1
        clone_cn_profiles = Matrix(clones_data[:, 2:end])
        ctdna = Vector{Float64}(ctdna_data[:, 1])
        scale = 1.0  

        log_potential = CtDNALogPotential(ctdna, clone_cn_profiles, num_clones, n, scale)
        reference_potential = default_reference(log_potential)

        time_taken = @elapsed begin
            pt = pigeons(
                target = log_potential,
                reference = reference_potential,
                record = [traces; record_default()],
                n_rounds = 4
            )
         report(pt)    
        end

        println("run complete for $ctdna_path. time taken: $time_taken seconds.")
    end
end

# ctdna_paths = ["data/ctdna.tsv","data/ctdna-500.tsv", "data/ctdna-1000.tsv","data/ctdna-5000.tsv"]
# clones_paths = ["data/2-clones-simple.tsv","data/2-clones-500.tsv", "data/2-clones-1000.tsv","data/2-clones-5000.tsv"]
ctdna_paths = ["data/ctdna.tsv"]
clones_paths = ["data/2-clones-simple.tsv"]
times = main(ctdna_paths, clones_paths)