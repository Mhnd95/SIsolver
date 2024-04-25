module SIsolver

# to test if the module is working and loaded properly
function hello()
    println("Hello, World!")
    x = 1
    y = 2 
    return x, y
end

using Revise
using CSV
using DataFrames
using Tables
# using Optimisers
using Optim
using LaTeXStrings


function read_data_file(filename::String)
    CSV.read(filename, DataFrame; header=false, skipto=6) |> Tables.matrix
end


function shape_factor()
    # F_s = 1 / V_ma * sum ( Area of ith surface open to flow / Distance of the ith surface from center of matrix ) from 1 to number of surfaces

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

function t_dimensionless(filename::String, θ::Float64)
    ExprmntData = read_data_file(filename)
    t_exprmnt = ExprmntData[:, 1]
    F_s = shape_factor()

    # dimensionless time equation
    # t_d = [ ( σ * cos(θ) F_s / µ_w ) * sqrt(k/Φ) ] t_exprmnt
    # σ = interfacial tension, θ = contact angle, µ_w = water viscosity, k = permeability, Φ = porosity
    σ = 14.12       # interfacial tension
    µ_w = 0.983     # water viscosity

    # k = ( 276.2 + 251 + 249.3 + 273.6 + 264 + 274.4 ) / 6            # permeability
    k = 276.2
    # Φ = ( 21.7 + 21.3 + 21.6 + 21.3 + 21.3 + 21.1 ) / 6            # porosity
    Φ = 21.7

    t_d = ( (σ * cos(θ) * F_s / µ_w) * sqrt(k / Φ) ) .* t_exprmnt
    

    return t_d
end

function R_experimental(filename::String)
    ExprmntData = read_data_file(filename)
    R_exprmnt = ExprmntData[:, 2]
    R_exprmnt_normalized = R_exprmnt ./ R_exprmnt[end]

    return R_exprmnt_normalized
end

function R_calculated(filename::String, θ::Float64, a::Vector{Float64}, λ::Vector{Float64})
    if length(a) != length(λ)
        throw(ArgumentError("Length of 'a' and 'λ' vectors must be equal."))
    end

    # Normalized calculated recovery factor equation:
    # (Rf/R_inf)_c = ( 1 -a[1]e^(-λ[1]*t_d)-a[2]e^(-λ[2]*t_d)-a[3]e^(-λ[3]*t_d))
    t_d = t_dimensionless(filename, θ)    

    R_calculated_normalized = 1 .- a[1] .* exp.(-λ[1] .* t_d) .- a[2] .* exp.(-λ[2] .* t_d) .- a[3] .* exp.(-λ[3] .* t_d)
    return R_calculated_normalized
end

function global_objective(params::Vector{Float64}, filename::String)
    a, λ, θ = params[1:3], params[4:6], params[7]
    R_exprmnt_normalized = R_experimental(filename)
    R_calculated_normalized = R_calculated(filename, θ, a, λ)
    
    if length(R_exprmnt_normalized) != length(R_calculated_normalized)
        throw(ArgumentError("Length of R_exprmnt_normalized and R_calculated_normalized vectors must be equal."))
    end

    # Residual sum of squares
    S_ri = sum((R_exprmnt_normalized .- R_calculated_normalized) .^ 2)
    return S_ri
end

# Optimization function for each file (Inner Optimization)
function optimize_theta(filename::String, a::Vector{Float64}, λ::Vector{Float64})
    obj_func = θ -> global_objective([a..., λ..., θ], filename)
    initial_theta = 0.0  # Initial guess for θ
    result = Optim.optimize(obj_func, initial_theta, NelderMead(), Optim.Options(show_trace=true))
    optimized_theta = result.minimizer
    return optimized_theta
end

# Outer optimization to adjust a and λ across multiple files
function optimize_params_across_files(filenames::Vector{String}, initial_params::Vector{Float64})
    obj_func = params -> sum(global_objective([params..., optimize_theta(filenames[i], params[1:3], params[4:6])], filenames[i]) for i in 1:length(filenames))
    result = Optim.optimize(obj_func, initial_params, NelderMead(), Optim.Options(show_trace=true))
    optimized_params = result.minimizer
    return optimized_params
end

function optimize_parameters(filename::String, initial_params::Vector{Float64})
    obj_func = params -> global_objective(params, filename)
    result = Optim.optimize(obj_func, initial_params, NelderMead(), Optim.Options(show_trace=true))
    optimized_params = result.minimizer

    output = "a = " * string(optimized_params[1:3]) * ", λ = " * string(optimized_params[4:6]) * ", θ = " * string(optimized_params[7])
    return output
end

#function optimize_parameters(filename::String, initial_params::Vector{Float64})
#    result = optimize(params -> global_objective(params, filename), initial_params)
#    optimized_params = result.minimizer

#    return optimized_params
#end

end