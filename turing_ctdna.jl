using Turing, Distributions, CSV, DataFrames, Pigeons
using InferenceReport

@model function ctDNA_model(ctdna, clone_cn_profiles, num_clones, n, scale)

    rho ~ Dirichlet(ones(num_clones))
    mu = Vector{Float64}(undef, n)
    total_sum = Vector{Float64}(undef, n)
    clone_cn_profiles = transpose(clone_cn_profiles)
    for i in 1:n
        total_sum[i] = sum(clone_cn_profiles[i, :] .* rho)
    end

    mean_total_sum = mean(total_sum)

    for i in 1:n
        mu[i] = log(total_sum[i]) - log(mean_total_sum)
    end

    degrees_of_freedom = 2
    for i in 1:n
        # no 3-parameter studentT distribution so we have to do scale and shift 
        ctdna[i] ~ TDist(degrees_of_freedom) * scale + mu[i]
    end
end


# load data
function load_data(ctdna_path, clones_path)
    ctdna_data = CSV.read(ctdna_path, DataFrame, delim='\t', header=false)
    clones_data = CSV.read(clones_path, DataFrame, delim='\t')

    return ctdna_data, clones_data
end

function main()
    ctdna_path = "data/ctdna.tsv"
    clones_path = "data/2-clones-simple.tsv"
    ctdna_data, clones_data = load_data(ctdna_path, clones_path)

    n = size(clones_data, 1)
    num_clones = size(clones_data, 2) - 1
    clone_cn_profiles = transpose(Matrix(clones_data[:, 2:end]))
    ctdna = Vector(ctdna_data[:, 1])
    scale = 1.0  # simplified scale for now

    # instantiate and run the model
    model = ctDNA_model(ctdna, clone_cn_profiles, num_clones, n, scale)
    turing_log_potential = TuringLogPotential(model)
    pt = pigeons(target = turing_log_potential,record = [traces; record_default()])
    println("Model run complete.")
    report(pt)
end


main()

