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

ssfm_result = SSFM(ψ₀, config)
#spectral_result = PseudoSpectral(ψ₀, config; solver = Tsit5())

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
