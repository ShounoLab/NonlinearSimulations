using FFTW
include("NLSE_utils.jl")


function Vₓ_NLSE(ψₜ :: Union{Vector{<:Real}, Vector{ComplexF64}})
    """
    potential for nonlinear Schrödinger equation
    """
    return -abs2.(ψₜ)
end


function Vₓ_barrier(
    inds :: Union{UnitRange{Int}, Vector{Int}},
    height :: Real,
    config :: NLSESettings)
    """
    barrer-type potentital
    """
    Vₓ = zeros(config.Ngrids)
    Vₓ[inds] .= height
    return Vₓ
end


function SSFM(
    ψ₀ :: Union{Vector{<:Real}, Vector{ComplexF64}},
    config :: NLSESettings)
    """
    Simulation with Split Step Fourier Method
    """

    ψ = zeros(ComplexF64, (config.Ngrids, config.Nsteps))
    Vₓ = zeros(ComplexF64, (config.Ngrids, config.Nsteps))
    ψ[:, 1] .= ψ₀
    Vₓ[:, 1] .= Vₓ_NLSE(ψ[:, 1])
    #Vₓ[:, 1] .= Vₓ_barrier(200:250, 5e10, config)
    for t in 2:config.Nsteps
        # potential valuse at time t
        #Vₓ[:, t] .= Vₓ_barrier(200:250, 5e10, config)
        Vₓ[:, t] .= Vₓ_NLSE(ψ[:, t - 1])
        # Solve the ODE about potential
        @views ψ₁ = exp.(im * config.Δt / 2 .* Vₓ[:, t]) .* ψ[:, t - 1]

        # Map to Fourier space
        Fψ₁ = fftshift(fft(ψ₁))

        # Time evolution in Fourier space
        Fψ₂ = exp.(im * config.Δt / 2 .* config.wavenums .^ 2) .* Fψ₁

        # Return to original space
        ψ[:, t] = ifft(fftshift(Fψ₂))
    end
    return NLSEResult(ψ, Vₓ, config)
end
