using Plots

include("ssfm.jl")
include("pseudospectral.jl")


### NLSE simulation by SSFM
Ngrids = 512 # Number of Fourier modes
L = 30.0 # Space period
Δt = 2π / 21
t_end = 2π
Nsteps = round(Int, t_end / Δt)

config = NLSESettings(Nsteps, Δt, t_end, Ngrids, L)

# initial state of the wave function
ψ₀ = 2.0 * sech.(config.gridpoints)

#result = SSFM(ψ₀, config)
result = PseudoSpectral(ψ₀, config)


# output a GIF animation
ymax = maximum(hcat(abs.(result.ψ), result.Vₓ))
anim = @animate for t in 1:result.config.Nsteps
    p = plot(result.config.gridpoints,
             abs.(result.ψ[:, t]),
             ylims = (0, ymax),
             label = "|psi|", lw = 2)
    plot!(result.config.gridpoints,
          result.Vₓ[:, t],
          label = "potential", lw = 2)
end

gif(anim, "nlse_ssfm.gif", fps = 60)

