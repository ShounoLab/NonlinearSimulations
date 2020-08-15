using DifferentialEquations
using FFTW
include("NLSE_utils.jl")


function spectral_evolution(Fψ, p, t)
    """
    time evolution function of NLSE on the Fourier space
    """
    ψ = ifft(Fψ)
    return - (im / 2) .* (fftshift(p.wavenums) .^ 2) .* Fψ .+ im .* fft(abs2.(ψ) .* ψ)
end


function PseudoSpectral(
    ψ₀ :: Union{Vector{<:Real}, Vector{ComplexF64}},
    config :: NLSESettings;
    solver = Tsit5())
    """
    Simulation with Pseudo-spectral Method
    """

    Fψ₀ = fft(ψ₀)
    prob = ODEProblem(spectral_evolution, Fψ₀, (0, config.t_end), config)
    sol = solve(prob, solver, dt = Δt, adaptive = true, reltol = 1e-8, abstol = 1e-8)
    Fψ = sol(0:config.Δt:config.t_end)

    ψ = Matrix{ComplexF64}(undef, (config.Ngrids, config.Nsteps))
    Vₓ = Matrix{ComplexF64}(undef, (config.Ngrids, config.Nsteps))
    for t in 1:config.Nsteps
        ψ[:, t] = ifft(Fψ[:, t])
        Vₓ[:, t] = abs2.(ψ[:, t])
    end

    return NLSEResult(ψ, Vₓ, config)
end
