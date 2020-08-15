using SparseArrays
include("Burgers_utils.jl")


function BurgersCrankNicolson_step(
    ψₜ :: Vector{<:Real},
    config :: BurgersSettings)
    """
    Unit-time evolution with Crank-Nicolson Method
    """

    p = config.Δt / (4 * config.Δx)
    r = ν * config.Δt / (2 * config.Δx ^ 2)

    aₜ = [- (r + p * ψₜ[n]) for n in 1:config.Ngrids]
    cₜ = [- r + p * ψₜ[n] for n in 1:config.Ngrids]
    bₜ = Vector{Float64}(undef, config.Ngrids)

    # periodic boundary conditions
    bₜ[1] = 2 * r + 1 + p * ψₜ[2] - p * ψₜ[config.Ngrids]
    bₜ[config.Ngrids] = 2 * r + 1 + p * ψₜ[1] - p * ψₜ[config.Ngrids - 1]
    for n in 2:(config.Ngrids - 1)
        bₜ[n] = 2 * r + 1 + p * ψₜ[n + 1] - p * ψₜ[n - 1]
    end

    Aₜ = spzeros(Float64, config.Ngrids, config.Ngrids)
    Bₜ = spzeros(Float64, config.Ngrids, config.Ngrids)
    for n in 1:config.Ngrids
        n₋ = ifelse(n == 1, config.Ngrids, n - 1)
        n₊ = ifelse(n == config.Ngrids, 1, n + 1)
        Aₜ[n, n₋] = aₜ[n]
        Aₜ[n, n] = bₜ[n]
        Aₜ[n, n₊] = cₜ[n]

        Bₜ[n, n₋] = r
        Bₜ[n, n] = 1 - 2 * r
        Bₜ[n, n₊] = r
    end

    return Aₜ \ (Bₜ * ψₜ)
end

function BurgersCrankNicolson(
    ψ₀ :: Vector{<:Real},
    config :: BurgersSettings)
    """
    Simulation with Crank-Nicolson Method
    """

    Ψ = Matrix{Float64}(undef, (config.Ngrids, config.Nsteps + 1))
    Ψ[:, 1] .= ψ₀
    for t in 1:config.Nsteps
        Ψ[:, t + 1] = BurgersCrankNicolson_step(Ψ[:, t], config)
    end

    return BurgersResult(Ψ, config)
end
