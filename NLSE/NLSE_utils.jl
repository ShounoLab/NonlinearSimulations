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


struct NLSEResult
    """
    Results of the NLSE simulation
    ψ (Matrix{Complex{Float64}}) : wave functions
    Vₓ (Array{Real}) : potentials
    config (NLSESettings) : configurations of NLSE
    """
    ψ :: Matrix{ComplexF64}
    Vₓ :: Array{Real}
    config :: NLSESettings
end
