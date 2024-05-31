using Pigeons
using Random
using StatsFuns
using CUDA

struct MyLogPotential
    n_trials::Int
    n_successes::Int
end

function (log_potential::MyLogPotential)(x)
    p1, p2 = x
    if !(0 < p1 < 1) || !(0 < p2 < 1)
        return -Inf64
    end
    p = p1 * p2
    return StatsFuns.binomlogpdf(log_potential.n_trials, p, log_potential.n_successes)
end

# e.g.:
my_log_potential = MyLogPotential(100, 50)
my_log_potential([0.5, 0.5])

# Pigeons.initialization(::MyLogPotential, ::AbstractRNG, ::Int) = [0.5, 0.5]
Pigeons.initialization(::MyLogPotential, ::AbstractRNG, ::Int) = [0.5, 0.5]


pt = pigeons(
        target = MyLogPotential(100, 50),
        reference = MyLogPotential(100, 100),
        record = [traces; record_default()]
    )

report(pt)