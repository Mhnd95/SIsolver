module SIsolver

function hello()
    println("Hello, World!")
end

using CSV
using DataFrames
using Tables
using Optimisers
using LaTeXStrings


function read_data_file(filename::String)
    CSV.read(filename; header=false, skipto=6, DataFrame) |> Tables.matrix
end


function shape_factor(volume::Float64, areas::Vector{Float64}, distances::Vector{Float64})
    # F_s = 1 / V_ma * sum ( Area of ith surface open to flow / Distance of the ith surface from center of matrix ) from 1 to number of surfaces

    if length(areas) != length(distances)
        throw(ArgumentError("Length of areas and distances vectors must be equal."))
    end
    areas = 2.0
    s = length(areas)
    sum_A_over_d = 0.0
    
    # Sample volumetric data, needs to be realistic core plug scale 
    volume = 0.5 * 1.5 ^ 2 * pi * 2.8
    areas = [0.5 * 1.5 ^ 2 * pi, 2.8 * 2 * pi * 1.5]
    distances = [1.5 , 1.4]

    for i in 1:s
        sum_A_over_d += areas[i] / distances[i]
    end

    F_s = (1 / volume) * sum_A_over_d
    return F_s
end

function experimental_results(filename::String, θ::Float64, F_s::Float64, σ::Float64, µ_w::Float64, k::Float64, Φ::Float64)
    ExprmntData = read_data_file(filename)
    t_exprmnt = ExprmntData[:, 1]
    R_exprmnt = ExprmntData[:, 2]
    
    # dimensionless time equation
    # t_D = [ ( σ * cos(θ) F_s / µ_w ) * sqrt(k/Φ) ] t_exprmnt
    # σ = interfacial tension, θ = contact angle, µ_w = water viscosity, k = permeability, Φ = porosity
    σ = 14.12       # interfacial tension
    µ_w = 0.983     # water viscosity

    # k = ( 276.2 + 251 + 249.3 + 273.6 + 264 + 274.4 ) / 6            # permeability
    k = 276.2
    # Φ = ( 21.7 + 21.3 + 21.6 + 21.3 + 21.3 + 21.1 ) / 6            # porosity
    Φ = 21.7

    t_dimensionless = ( (σ * cos(θ) * F_s / µ_w) * sqrt(k / Φ) ) .* t_exprmnt
    R_exprmnt_normalized = R_exprmnt ./ R_exprmnt[end]

    return t_dimensionless, R_exprmnt_normalized
end

function R_calculated(t_dimensionless::Vector{Float64}, a::Vector{Float64}, λ::Vector{Float64})
    # Normalized calculated recovery factor equation:
    # (Rf/R_inf)_c = ( 1 -a[1]e^(-λ[1]*t_dimensionless)-a[2]e^(-λ[2]*t_dimensionless)-a[3]e^(-λ[3]*t_dimensionless))

    if length(a) != length(λ)
        throw(ArgumentError("Length of 'a' and 'λ' vectors must be equal."))
    end

    R_calculated_normalized = 1 .- a[1] .* exp.(-λ[1] .* t_dimensionless) .- a[2] .* exp.(-λ[2] .* t_dimensionless) .- a[3] .* exp.(-λ[3] .* t_dimensionless)
    return R_calculated_normalized
end

function global_objective(R_exprmnt_normalized::Vector{Float64}, R_calculated_normalized::Vector{Float64})
    if length(R_exprmnt_normalized) != length(R_calculated_normalized)
        throw(ArgumentError("Length of R_exprmnt_normalized and R_calculated_normalized vectors must be equal."))
    end

    # Global objective function that will be used to calculate the error
    S_ri = sum((R_exprmnt_normalized .- R_calculated_normalized) .^ 2)
    # How do I make sure that point 1 in R_exprmnt_normalized is compared to R_calculated_normalized?
    # correlated_data = hcat(t_dimensionless, R_exprmnt_normalized, R_calculated_normalized)
    O_g = S_ri
    return O_g
end

function optimization_objective(params::Vector{Float64}, filename::String, F_s::Float64, σ::Float64, µ_w::Float64, k::Float64, Φ::Float64)
    a = params[1:3]
    λ = params[4:6]
    θ = params[7]

    t_dimensionless, R_exprmnt_normalized = experimental_results(filename, θ, F_s, σ, µ_w, k, Φ)
    R_calculated_normalized = R_calculated(t_dimensionless, a, λ)
    O_g = global_objective(R_exprmnt_normalized, R_calculated_normalized)

    return O_g
end

function optimize_parameters(filename::String, F_s::Float64, σ::Float64, µ_w::Float64, k::Float64, Φ::Float64, initial_params::Vector{Float64})
    result = optimize(params -> optimization_objective(params, filename, F_s, σ, µ_w, k, Φ), initial_params)
    optimized_params = result.minimizer

    return optimized_params
end

end