struct BurgersSettings
    """
    Settings for simulating Burgers Equation.
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
    ν :: Real
end

function BurgersSettings(
    Nsteps :: Int64,
    Δt :: Real,
    t_end :: Real,
    Ngrids :: Int64,
    Xsize :: Real,
    ν :: Real
    )
    """
    outer constructor of BurgersSettings
    """
    Δx = Xsize / Ngrids
    space_inds = (-Ngrids / 2):(Ngrids / 2 - 1)
    xs = space_inds .* Δx
    ks = 2 .* space_inds .* π / L
    return BurgersSettings(Nsteps, Δt, t_end, Ngrids, Xsize, Δx, collect(xs), collect(ks), ν)
end


struct BurgersResult
    """
    Results of the Burgers simulation
    ψ (Matrix{Float64}) : wave functions
    config (BurgersSettings) : configurations of Burgers Equation
    """
    ψ :: Matrix{Float64}
    config :: BurgersSettings
end
