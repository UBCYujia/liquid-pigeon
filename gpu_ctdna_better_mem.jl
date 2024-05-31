using CUDA
using Pigeons
using Random
using Statistics
using Distributions
using CSV
using DataFrames
using InferenceReport
using SpecialFunctions

CUDA.allowscalar(false)  


#global gpu data
global GLOBAL_CTDNA = Ref{CuArray{Float32, 1}}()
global GLOBAL_CLONE_CN_PROFILES = Ref{CuArray{Float32, 2}}()
global GLOBAL_RHO = CuArray{Float32, 1}(undef, 2)

function load_data_to_gpu(ctdna_path, clones_path)
    start_time = time_ns()

    ctdna_data = CSV.read(ctdna_path, DataFrame, delim='\t', header=false, types=[Float32])
    clones_data = CSV.read(clones_path, DataFrame, delim='\t')

    GLOBAL_CTDNA[] = CuArray(Float32.(Vector(ctdna_data[:, 1])))
    GLOBAL_CLONE_CN_PROFILES[] = CuArray(Float32.(Matrix(clones_data[:, 2:end])))

    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 1e9
    #println("Execution of load_data_to_gpu took $elapsed_time seconds")
end


struct CtDNALogPotential
    num_clones::Int
    n::Int
    scale::Float32
end

function log_t_pdf(x, v)
    result = - ((v + 1) / 2) .* log.(1 .+ (x .^ 2) ./ v)
    return result
end 

# function (log_potential::CtDNALogPotential)(params)
#     rho = CuArray(Float32.(params)) 
#     if any(rho .< 0 .|| rho .> 1) || abs(sum(rho) - 1) > 1e-5
#         return -Inf
#     end

#     total_sum = GLOBAL_CLONE_CN_PROFILES[] * rho
#     mean_total_sum = mean(CUDA.reduce(+, total_sum) / length(total_sum))

#     mu = log.(total_sum) .- log(mean_total_sum)
#     degrees_of_freedom = 2
#     scaled_mu = mu * log_potential.scale
#     log_likelihoods = log_t_pdf((GLOBAL_CTDNA[] .- scaled_mu) / log_potential.scale, degrees_of_freedom)
#     log_likelihood = CUDA.reduce(+, log_likelihoods)
    
#     return log_likelihood
# end

function (log_potential::CtDNALogPotential)(params)
    

    
    start_time = time_ns()
    if any(x -> x < 0 || x > 1, params) || abs(sum(params)  - 1) > 1e-5
        # if any(x -> x < 0 || x > 1, rho) || sum(rho) != 1
            return -Inf  #ensure rho is valid
    end

    
    copyto!(GLOBAL_RHO, Float32.(params)) 
    total_sum = GLOBAL_CLONE_CN_PROFILES[] * GLOBAL_RHO

    
    mean_total_sum = mean(CUDA.reduce(+, total_sum) / length(total_sum))
    
    mu = log.(total_sum) .- log(mean_total_sum)
    end_time = time_ns()
    
    degrees_of_freedom = 2
    scaled_mu = mu * log_potential.scale
    log_likelihoods = log_t_pdf((GLOBAL_CTDNA[] .- scaled_mu) / log_potential.scale, degrees_of_freedom)
    log_likelihood = CUDA.reduce(+, log_likelihoods)

    
    elapsed_time = (end_time - start_time) / 1e9
    #println("Execution of log_potential took $elapsed_time seconds")
    
    return log_likelihood
end


function Pigeons.initialization(log_potential::CtDNALogPotential, rng::AbstractRNG, ::Int)

    start_time = time_ns()

    alpha = 1.0  
    rho = rand(rng, Dirichlet(log_potential.num_clones, alpha))  

    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 1e9
    #println("Execution of pigeons_init took $elapsed_time seconds")
    return rho  # cannot convert to cuarray
end

function Pigeons.sample_iid!(log_potential::CtDNALogPotential, replica, shared)
    start_time = time_ns()

    rng = replica.rng
    new_state = rand(rng, Dirichlet(log_potential.num_clones, 1.0)) #to cuda rand?

    @assert abs(sum(new_state) - 1) < 1e-5 "density not 1!"

    replica.state = new_state

    end_time = time_ns()
    elapsed_time = (end_time - start_time) / 1e9
    #println("Execution of pigeons_sampleiid took $elapsed_time seconds")
end

function default_reference(log_potential::CtDNALogPotential)
    return CtDNALogPotential(log_potential.num_clones, log_potential.n, 1)
end

function main(ctdna_paths, clones_paths)
    times = Float64[]
    load_data_to_gpu(ctdna_paths[1], clones_paths[1]) 

    for (ctdna_path, clones_path) in zip(ctdna_paths, clones_paths)
        println("processing: $ctdna_path and $clones_path")
        n = length(GLOBAL_CTDNA[])
        num_clones = size(GLOBAL_CLONE_CN_PROFILES[], 2)
        scale = 1.0

        log_potential = CtDNALogPotential(num_clones, n, scale)
        reference_potential = default_reference(log_potential)  

        profile_res = @elapsed begin
            pt = pigeons(
                target = log_potential,
                reference = reference_potential,
                record = [traces; record_default()],
                n_rounds = 10,
                n_chains = 10
            )
        end
        #push!(times, time_taken)
        println(profile_res)
    end
    return times
end

ctdna_paths = ["data/ctdna-5000.tsv"]
clones_paths = ["data/2-clones-5000.tsv"]
times = main(ctdna_paths, clones_paths)