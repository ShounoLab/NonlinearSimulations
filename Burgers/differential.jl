include("./Burgers_utils.jl")


function Differential(
    ψ₀ :: Vector{<:Real},
    config :: BurgersSettings;
    Δt_intrinsic = 1e-4)
    """
    Simulation with Differential Method
    """

    t_ary_intrinsic = collect(0:Δt_intrinsic:config.t_end)
    t_ary = collect(0:config.Δt:config.t_end)
    Nsteps_intrinsic = length(t_ary_intrinsic)
    t_ary_idx = [findmin(abs.(t_ary[i] .- t_ary_intrinsic))[2] for i in 1:length(t_ary)]
    Ψ = Matrix{Float64}(undef, (config.Ngrids, Nsteps_intrinsic))
    Ψ[:, 1] .= ψ₀
    for t in 2:Nsteps_intrinsic
        for n in 1:config.Ngrids
            if n == 1
                Ψ[n, t] = burgers_rhs(Ψ[end, t - 1], Ψ[n, t - 1], Ψ[n + 1, t - 1],
                                      Δt_intrinsic, config.Δx, config.ν)
            elseif n == config.Ngrids
                Ψ[n, t] = burgers_rhs(Ψ[n - 1, t - 1], Ψ[n, t - 1], Ψ[1, t - 1],
                                      Δt_intrinsic, config.Δx, config.ν)
            else
                Ψ[n, t] = burgers_rhs(Ψ[n - 1, t - 1], Ψ[n, t - 1], Ψ[n + 1, t - 1],
                                      Δt_intrinsic, config.Δx, config.ν)
            end
        end
    end

    return BurgersResult(Ψ[:, t_ary_idx], config)
end


