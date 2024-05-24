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
    ctdna_data = CSV.read(ctdna_path, DataFrame, delim='\t', header=false)
    clones_data = CSV.read(clones_path, DataFrame, delim='\t')
    return ctdna_data, clones_data
end

function default_reference(log_potential::CtDNALogPotential)
    neutral_ctdna = ones(log_potential.n)
    return CtDNALogPotential(neutral_ctdna, log_potential.clone_cn_profiles, log_potential.num_clones, log_potential.n, log_potential.scale)
end

function main()
    ctdna_path = "data/ctdna.tsv"
    clones_path = "data/2-clones-simple.tsv"
    ctdna_data, clones_data = load_data(ctdna_path, clones_path)

    n = size(clones_data, 1)
    num_clones = size(clones_data, 2) - 1
    clone_cn_profiles = Matrix(clones_data[:, 2:end])
    ctdna = Vector{Float64}(ctdna_data[:, 1])
    scale = 1.0  #simplified scale

    log_potential = CtDNALogPotential(ctdna, clone_cn_profiles, num_clones, n, scale)
    reference_potential = default_reference(log_potential)

    pt = pigeons(
        target = log_potential,
        reference = reference_potential,
        record = [traces; record_default()]
    )
    report(pt)

    println("Model run complete.")
end

main()
