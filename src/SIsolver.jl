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
using Plots


function read_data_file(filename::String)
    CSV.read(filename, DataFrame; header=false, skipto=6) |> Tables.matrix
end


function shape_factor()
    # F_s = 1 / V_ma * sum ( Area of ith surface open to flow / Distance of the ith surface from center of matrix ) from 1 to number of surfaces
    areas = [0.5 * 1.5^2 * pi, 2.8 * 2 * pi * 1.5]
    distances = [1.5, 1.4]
    volume = 0.5 * 1.5^2 * pi * 2.8
    sum_A_over_d = sum(areas[i] / distances[i] for i in 1:length(areas))
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

    R_calculated_normalized = 1 .- sum(a[i] .* exp.(-λ[i] .* t_d) for i in 1:length(a))
    return t_d, R_calculated_normalized
end

function global_objective(params::Vector{Float64}, filename::String)
    a, λ, θ = params[1:3], params[4:6], params[7]
    R_exprmnt_normalized = R_experimental(filename)
    _, R_calculated_normalized = R_calculated(filename, θ, a, λ)

    if length(R_exprmnt_normalized) != length(R_calculated_normalized)
        println("R_exprmnt_normalized length: ", length(R_exprmnt_normalized))
        println("R_calculated_normalized length: ", length(R_calculated_normalized))
        throw(ArgumentError("Length of R_exprmnt_normalized and R_calculated_normalized vectors must be equal."))
    end

    # Residual sum of squares
    S_ri = sum((R_exprmnt_normalized .- R_calculated_normalized) .^ 2)
    return S_ri
end

# Optimization function for each file (Inner Optimization)
function optimize_theta(filename::String, a::Vector{Float64}, λ::Vector{Float64})
    obj_func = θ -> global_objective(vcat(a, λ, θ), filename)
    initial_theta = [0.0]  # Initial guess for θ
    lower_bound = [0.0]
    upper_bound = [π]
    result = Optim.optimize(obj_func, lower_bound, upper_bound, initial_theta, Fminbox(LBFGS()), Optim.Options(show_trace=true))
    optimized_theta = result.minimizer[1]
    return optimized_theta
end

# Outer optimization to adjust a and λ across multiple files
function optimize_params_across_files(filenames::Vector{String})
    initial_params = [0.33, 0.33, 0.33, 0.01, 0.01, 0.01]
    obj_func = params -> begin
        a, λ = params[1:3], params[4:6]
        total_error = 0.0
        for i in 1:length(filenames)
            θ = optimize_theta(filenames[i], a, λ)
            total_error += global_objective(vcat(a, λ, θ), filenames[i])
        end
        return total_error
    end

    lower_bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Lower bounds for a and λ
    upper_bounds = [Inf, Inf, Inf, Inf, Inf, Inf] # No upper bounds for a and λ
    
    result = Optim.optimize(obj_func, lower_bounds, upper_bounds, initial_params, Fminbox(LBFGS()), Optim.Options(show_trace=true,))
    optimized_params = result.minimizer
    return optimized_params
end

function plot_results(file_pattern::String, save_path::String)
    # Find all relevant files in the folder
    data_folder = "data"
    filenames = filter(x -> occursin(file_pattern, x), readdir(data_folder))
    filenames = [joinpath(data_folder, x) for x in filenames]
    
    output_folder = "output"
    save_path = joinpath(output_folder, save_path)

    optimized_params = optimize_params_across_files(filenames)
    a, λ = optimized_params[1:3], optimized_params[4:6]

    p_error = plot(title="Error vs Time for Multiple Files", xlabel="Time", ylabel="Error")
    p_R_vs_td = plot(title="R_calculated vs t_d for Multiple Files", xlabel="t_d", ylabel="R_calculated")

    θ_values = Float64[]
    for filename in filenames
        θ = optimize_theta(filename, a, λ)
        push!(θ_values, θ)
        ExprmntData = read_data_file(filename)
        t_exprmnt = ExprmntData[:, 1]
        R_exprmnt_normalized = R_experimental(filename)
        t_d, R_calculated_normalized = R_calculated(filename, θ, a, λ)
        
        if length(R_exprmnt_normalized) != length(R_calculated_normalized)
            println("Filename: ", filename)
            println("R_exprmnt_normalized length: ", length(R_exprmnt_normalized))
            println("R_calculated_normalized length: ", length(R_calculated_normalized))
            throw(ArgumentError("Length of R_exprmnt_normalized and R_calculated_normalized vectors must be equal."))
        end

        error = R_exprmnt_normalized .- R_calculated_normalized
        
        plot!(p_error, t_exprmnt, error, label=filename, lw=2)
        plot!(p_R_vs_td, t_d, R_calculated_normalized, label=filename, lw=2)
    end

    savefig(p_error, save_path * "_error.png")
    savefig(p_R_vs_td, save_path * "_R_vs_td.png")
    
    # Convert θ values to degrees
    θ_values_deg = rad2deg.(θ_values)
    
    # Create a DataFrame to display optimized parameters
    df = DataFrame(
        File = filenames,
        a1 = a[1],
        a2 = a[2],
        a3 = a[3],
        λ1 = λ[1],
        λ2 = λ[2],
        λ3 = λ[3],
        θ_deg = θ_values_deg
    )
    
    println("Optimized Parameters:")
    println(df)
end


function optimize_parameters(filename::String)
    initial_params = [0.33, 0.33, 0.33, 0.01, 0.01, 0.01, 0.0]
    obj_func = params -> global_objective(params, filename)
    #lower_bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # Lower bounds for a and λ
    #upper_bounds = [Inf, Inf, Inf, Inf, Inf, Inf] # No upper bounds for a and λ
    
    #result = Optim.optimize(obj_func, lower_bounds, upper_bounds, initial_params, Fminbox(NelderMead()), Optim.Options(show_trace=true))
    result = Optim.optimize(obj_func, initial_params, NelderMead(), Optim.Options(show_trace=true))
    optimized_params = result.minimizer

    output = "a = " * string(optimized_params[1:3]) * ", λ = " * string(optimized_params[4:6]) * ", θ = " * string(optimized_params[7])
    return output
end

function plot_error_vs_time(filename::String, θ::Float64, a::Vector{Float64}, λ::Vector{Float64}, save_path::String)
    ExprmntData = read_data_file(filename)
    t_exprmnt = ExprmntData[:, 1]
    R_exprmnt_normalized = R_experimental(filename)
    R_calculated_normalized = R_calculated(filename, θ, a, λ)
    
    if length(R_exprmnt_normalized) != length(R_calculated_normalized)
        throw(ArgumentError("Length of R_exprmnt_normalized and R_calculated_normalized vectors must be equal."))
    end

    error = R_exprmnt_normalized .- R_calculated_normalized
    
    p = plot(t_exprmnt, error, xlabel="Time", ylabel="Error", title="Error vs Time", label="Error", lw=2)
    savefig(p, save_path)
end

end