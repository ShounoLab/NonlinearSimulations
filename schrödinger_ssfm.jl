using FFTW
using Plots


struct NLSESettings
    """
    Settings for simulating Nonlinear Schrödinger Equation.
    Nsteps (Int64) : number of timepoints
    Δt (Real) : periods of timesteps
    t_end (Real) : endpoint
    Ngrids (Int64) : number of grids
    Xsize (Real) : system size
    Δx (Real) : grid size
    gridpoints (Vector{Real}) : locations of grid points
    wavenums (Vector{Real}) : wave numbers for FFT """
    Nsteps :: Int
    Δt :: Real
    t_end :: Real
    Ngrids :: Int
    Xsize :: Real
    Δx :: Real
    gridpoints :: Vector{Real}
    wavenums :: Vector{Real}
end

function NLSESettings(
    Nsteps :: Int64,
    Δt :: Real,
    t_end :: Real,
    Ngrids :: Int64,
    Xsize :: Real
    )
    """
    outer constructor of NLSESettings
    """
    Δx = Xsize / Ngrids
    space_inds = (-Ngrids / 2):(Ngrids / 2 - 1)
    xs = space_inds .* Δx
    ks = 2 .* space_inds .* π / L
    return NLSESettings(Nsteps, Δt, t_end, Ngrids, Xsize, Δx, collect(xs), collect(ks))
end


struct SSFMResult
    """
    Results of the Split Step Fourier Method
    ψ (Matrix{Complex{Float64}}) : wave functions
    Vₓ (Array{Real}) : potentials
    config (NLSESettings) : configurations of NLSE
    """
    ψ :: Matrix{ComplexF64}
    Vₓ :: Array{Real}
    config :: NLSESettings
end



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
    return SSFMResult(ψ, Vₓ, config)
end


Ngrids = 512 # Number of Fourier modes
L = 30.0 # Space period
Δt = 2π / 21
t_end = 2π
Nsteps = round(Int, t_end / Δt)

config = NLSESettings(Nsteps, Δt, t_end, Ngrids, L)

# initial state of the wave function
ψ₀ = 2.0 * sech.(config.gridpoints)

ssfm_result = SSFM(ψ₀, config)

heatmap(abs.(ssfm_result.ψ))

anim = @animate for t in 1:Nsteps
    p = plot(ssfm_result.config.gridpoints,
             abs.(ssfm_result.ψ[:, t]),
             ylims = (0, maximum(abs.(ssfm_result.ψ))),
             label = :none, lw = 2)
    plot!(ssfm_result.config.gridpoints,
          ssfm_result.Vₓ[:, t] ./ maximum(ssfm_result.Vₓ[:, t]) .*
          maximum(abs.(ssfm_result.ψ)),
          label = "potential (scaled)", lw = 2)
end

gif(anim, "nlse.gif", fps = 30)
