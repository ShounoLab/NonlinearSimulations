using Plots

include("./Burgers/cranknicolson.jl")
#.nclude("./Burgers/differential.jl")


### Burgers' equation
Ngrids = 256 # Number of Fourier modes
L = 30.0 # Space period
Δt = 0.1
t_end = 30.0
Nsteps = round(Int, t_end / Δt)
ν = 0.1

config = BurgersSettings(Nsteps, Δt, t_end, Ngrids, L, ν)

# initial state of the wave function
ψ₀ = exp.(- (config.gridpoints .+ 2) .^ 2)

result = BurgersCrankNicolson(ψ₀, config)


# output a GIF animation
ymin, ymax = minimum(result.ψ), maximum(result.ψ)
anim = @animate for t in 1:result.config.Nsteps
    p = plot(result.config.gridpoints,
             abs.(result.ψ[:, t]),
             ylims = (ymin, ymax),
             label = "|psi|", lw = 2)
end

gif(anim, "burgers_cranknicolson.gif", fps = 60)
