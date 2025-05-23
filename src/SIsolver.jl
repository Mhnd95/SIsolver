module SIsolver

# to test if the module is loaded and working properly
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
using Optim
using LaTeXStrings
using Plots
using StatsPlots
using Dates
using Random
using Printf

function read_data_file(filename::String)
    CSV.read(filename, DataFrame; header=false, skipto=6) |> Tables.matrix
end

function shape_factor()
    areas = [1.5^2 * pi, 1.5^2 * pi, 2.8 * 2 * pi * 1.5]
    distances = [1.4, 1.4, 1.5]
    volume = 1.5^2 * pi * 2.8
    sum_A_over_d = sum(areas[i] / distances[i] for i in 1:length(areas))
    F_s = (1 / volume) * sum_A_over_d

    return F_s
end

function t_dimensionless(filename::String, θ::Float64)
    ExprmntData = read_data_file(filename)
    t_exprmnt = ExprmntData[:, 1]
    F_s = shape_factor()

    σ = 14.12       # interfacial tension
    µ_w = 0.983     # water viscosity
    k = 276.2       # permeability
    Φ = 21.7        # porosity

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
    t_d = t_dimensionless(filename, θ)    

    R_calculated_normalized = 1 .- sum(a[i] .* exp.(-λ[i] .* t_d) for i in 1:length(a))
    return t_d, R_calculated_normalized
end

function global_objective(params::Vector{Float64}, filenames::Vector{String})
    n_files = length(filenames)
    a = params[1:3]
    λ = params[4:6]
    θ_values = params[7:end]

    # Normalize a to ensure the sum is 1
    a = a / sum(a)

    total_error = 0.0
    for i in 1:n_files
        θ = θ_values[i]
        R_exprmnt_normalized = R_experimental(filenames[i])
        _, R_calculated_normalized = R_calculated(filenames[i], θ, a, λ)

        if length(R_exprmnt_normalized) != length(R_calculated_normalized)
            println("R_exprmnt_normalized length: ", length(R_exprmnt_normalized))
            println("R_calculated_normalized length: ", length(R_calculated_normalized))
            throw(ArgumentError("Length of R_exprmnt_normalized and R_calculated_normalized vectors must be equal."))
        end

        total_error += sum((R_exprmnt_normalized .- R_calculated_normalized) .^ 2)
    end

    return total_error
end

# Optimization function for each file (Inner Optimization)
function optimize_theta(filename::String, a::Vector{Float64}, λ::Vector{Float64})
    obj_func = θ -> global_objective(vcat(a, λ, θ), filename)
    initial_theta = [0.0]  # Initial guess for θ
    lower_bound = [0.0]
    upper_bound = [π]
    result = Optim.optimize(obj_func, lower_bound, upper_bound, initial_theta, Fminbox(NelderMead()), Optim.Options(show_trace=true))
    optimized_theta = result.minimizer[1]
    return optimized_theta
end

# Outer optimization to adjust a and λ across multiple files
function optimize_params_across_files(filenames::Vector{String}, max_iter::Int=1000)
    n_files = length(filenames)
    initial_params = [0.33, 0.33, 0.33, 0.01, 0.01, 0.01]
    initial_θ = fill(0.0, n_files)
    initial_params = vcat(initial_params, initial_θ)

    obj_func = params -> begin
        # Ensure the sum of the first three parameters (a) is 1
        a = params[1:3] / sum(params[1:3])
        remaining_params = params[4:end]
        normalized_params = vcat(a, remaining_params)
        global_objective(normalized_params, filenames)
    end

    lower_bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, fill(0.0, n_files)...]
    upper_bounds = [Inf, Inf, Inf, Inf, Inf, Inf, fill(π, n_files)...]
    
    options = Optim.Options(store_trace=true, show_trace=true, iterations=max_iter)
    
    result = Optim.optimize(obj_func, lower_bounds, upper_bounds, initial_params, Fminbox(NelderMead()), options)
    optimized_params = result.minimizer
    trace = result.trace

    # Ensure the sum of the first three parameters (a) is 1
    a = optimized_params[1:3] / sum(optimized_params[1:3])
    remaining_params = optimized_params[4:end]
    normalized_params = vcat(a, remaining_params)
    
    objective_history = [trace[i].value for i in 1:length(trace)]
    
    return normalized_params, objective_history
end

function save_plot_data(data, file_path::String)
    CSV.write(file_path, DataFrame(data))
end

function plot_results(file_pattern::String, save_path::String, max_iter::Int=1000)
    start_time = now()
    optimization_method = "Fminbox(NelderMead())"
    data_folder = "data"
    filenames = filter(x -> occursin(file_pattern, x), readdir(data_folder))
    filenames = [joinpath(data_folder, x) for x in filenames]
    
    output_folder = "output"
    save_path = joinpath(output_folder, save_path)

    # Colors for consistency across plots
    colors = distinguishable_colors(length(filenames))

    optimized_params, objective_history = optimize_params_across_files(filenames, max_iter)
    a, λ = optimized_params[1:3], optimized_params[4:6]
    θ_values = optimized_params[7:end]

    # Data storage
    plot_data = Dict{String, Any}()

    # R_vs_td plot data
    p_R_vs_td_data = DataFrame(Time_Dimensionless=Float64[], R_Calculated=Float64[], File=String[])
    p_R_vs_td = plot(
        title=L"R_{calculated} \, vs \, t_d \, for \, Multiple \, Files", 
        xlabel=L"Dimensionless \, Time \, (t_d)", 
        ylabel=L"R_{calculated}", 
        legend=:topright, 
        grid=false,
        background_color=:white,
        xlims=(0, 10000),
        ylims=(0, 1)
    )
    
    # Residuals plot data
    p_residuals_data = DataFrame(Time=Float64[], Residuals=Float64[], File=String[])
    p_residuals = plot(
        title=L"Residuals \, for \, Multiple \, Files", 
        xlabel=L"Time \, (t)", 
        ylabel=L"Residuals", 
        legend=:topright, 
        grid=false,
        background_color=:white,
        xlims=(0, 120),
        ylims=(-0.2, 0.2)
    )
    all_residuals = Dict{String, Vector{Float64}}()

    for (i, filename) in enumerate(filenames)
        θ = θ_values[i]
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
        residuals = R_exprmnt_normalized .- R_calculated_normalized
        all_residuals[filename] = residuals

        color = colors[i]

        plot!(p_R_vs_td, t_d, R_calculated_normalized, label=filename, lw=2, color=color)
        plot!(p_residuals, t_exprmnt, residuals, label=filename, lw=2, color=color)

        # Store R_vs_td data
        append!(p_R_vs_td_data, DataFrame(Time_Dimensionless=t_d, R_Calculated=R_calculated_normalized, File=repeat([filename], length(t_d))))

        # Store residuals data
        append!(p_residuals_data, DataFrame(Time=t_exprmnt, Residuals=residuals, File=repeat([filename], length(residuals))))

        # Individual model fit plot for each file
        clean_filename = replace(basename(filename), ".csv" => "")
        model_fit_title_string = L"Model \, Fit \, for \, " * LaTeXString(clean_filename)
        p_model_fit = plot(
            title=model_fit_title_string, 
            xlabel=L"Time \, (t)", 
            ylabel=L"R", 
            legend=:topright, 
            grid=false,
            background_color=:white,
            xlims=(0, 120),
            ylims=(0, 1)
        )
        plot!(p_model_fit, t_exprmnt, R_exprmnt_normalized, label=L"Experimental", lw=2, linestyle=:solid, color=color)
        plot!(p_model_fit, t_exprmnt, R_calculated_normalized, label=L"Calculated", lw=2, linestyle=:dash, color=color)
        savefig(p_model_fit, save_path * "_model_fit_" * replace(clean_filename, " " => "_") * ".svg")
    end

    # Save R_vs_td plot data
    save_plot_data(p_R_vs_td_data, save_path * "_R_vs_td_data.csv")
    
    # Save residuals plot data
    save_plot_data(p_residuals_data, save_path * "_residuals_data.csv")
    
    savefig(p_R_vs_td, save_path * "_R_vs_td.svg")
    savefig(p_residuals, save_path * "_residuals.svg")

    # Histogram of residuals
    p_hist_residuals = histogram(collect(Iterators.flatten(values(all_residuals))), bins=30, 
        title=L"Histogram \, of \, Residuals", 
        xlabel=L"Residuals", 
        ylabel=L"Frequency", 
        legend=false, 
        grid=false,
        background_color=:white
    )
    savefig(p_hist_residuals, save_path * "_hist_residuals.svg")
    
    # Clean filenames for labels
    clean_filename = filename -> replace(replace(basename(filename), r"\.csv" => ""), r"data/" => "")
    # Box plot of residuals
    residuals_per_file = [R_experimental(filename) .- R_calculated(filename, θ_values[i], a, λ)[2] for (i, filename) in enumerate(filenames)]
    residual_labels = [clean_filename(f) for f in filenames]

    p_box_residuals = boxplot(
        residuals_per_file,
        positions = 1:length(filenames),
        title=L"Box \, Plot \, of \, Residuals",
        ylabel=L"Residuals",
        xticks=(1:length(filenames), residual_labels),  # Correct x-axis labels
        legend=false,
        grid=false,
        background_color=:white,
        ylims=(-0.2, 0.2)
    )
    savefig(p_box_residuals, save_path * "_box_residuals.svg")

    # Save box plot data
    box_plot_data = DataFrame(File=String[], Residuals=Float64[])
    for (filename, residuals) in all_residuals
        append!(box_plot_data, DataFrame(File=repeat([replace(basename(filename), ".csv" => "")], length(residuals)), Residuals=residuals))
    end
    save_plot_data(box_plot_data, save_path * "_box_plot_residuals_data.csv")

    # Heatmap of parameter sensitivity
    sensitivity_data = [a λ]
    p_heatmap_sensitivity = heatmap(sensitivity_data, 
        title=L"Heatmap \, of \, Parameter \, Sensitivity", 
        xlabel=L"Parameters \, a", 
        ylabel=L"Parameters \, \lambda", 
        color=:viridis,
        background_color=:white
    )
    savefig(p_heatmap_sensitivity, save_path * "_heatmap_sensitivity.svg")
    
    θ_values_deg = rad2deg.(θ_values)
    end_time = now()
    execution_time = end_time - start_time
    total_seconds = Int(floor(Nanosecond(execution_time).value // 1_000_000_000))
    hours, rem = divrem(total_seconds, 3600)
    minutes, seconds = divrem(rem, 60)
    formatted_time = @sprintf("%02d:%02d:%02d", hours, minutes, seconds)
    df = DataFrame(
        File = filenames,
        a1 = a[1],
        a2 = a[2],
        a3 = a[3],
        λ1 = λ[1],
        λ2 = λ[2],
        λ3 = λ[3],
        θ_deg = θ_values_deg,
        OptimizationMethod = repeat([optimization_method], length(filenames)),
        ExecutionTime = repeat([formatted_time], length(filenames))
    )
    
    println("Optimized Parameters:")
    println(df)

    # Save the DataFrame to a CSV file
    results_file = save_path * "_optimized_parameters.csv"
    CSV.write(results_file, df)
    println("Optimized parameters saved to ", results_file)

    # Plot optimization objective function history
    p_objective_history = plot(
        1:length(objective_history), objective_history, 
        xlabel=L"Iteration", 
        ylabel=L"Objective \, Function \, Value", 
        title=L"Optimization \, Objective \, Function \, vs. \, Iteration", 
        lw=2, 
        label="Objective Function History",
        legend=:topright,
        grid=false,
        background_color=:white
    )
    savefig(p_objective_history, save_path * "_objective_history.svg")

    # Save optimization objective function history data
    save_plot_data(DataFrame(Iteration=1:length(objective_history), ObjectiveFunctionValue=objective_history), save_path * "_objective_history_data.csv")

    println("Execution Time: ", formatted_time)
end

# Optimization for one single file:
function optimize_parameters(filename::String)
    initial_params = [0.33, 0.33, 0.33, 0.01, 0.01, 0.01, 0.0]
    obj_func = params -> global_objective(params, filename)
    
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

function generate_synthetic_data(filename::String, params::Vector{Float64}, num_points::Int=100)
    # Known parameters
    a = params[1:3] / sum(params[1:3])  # Normalize a to sum to 1
    λ = params[4:6]
    θ = deg2rad(params[7])  # Convert degrees to radians
    
    # Experimental setup
    σ = 14.12             # Interfacial tension
    µ_w = 0.983           # Water viscosity
    k = 276.2             # Permeability
    Φ = 21.7              # Porosity
    F_s = shape_factor()  # Get shape factor from the same setup for simplicity

    # Generate synthetic time points
    t_exprmnt = range(1, stop=120, length=num_points)  # 1 to 120 minutes

    # Calculate dimensionless time
    t_d = ( (σ * cos(θ) * F_s / µ_w) * sqrt(k / Φ) ) .* t_exprmnt

    # Generate synthetic recovery factor data
    R_calculated = 1 .- sum(a[i] .* exp.(-λ[i] .* t_d) for i in 1:length(a))
    
    # Add noise to the data
    noise_level = 0.01  # 1% noise
    R_noisy = R_calculated .+ noise_level .* randn(length(R_calculated))
    
    # Ensure recovery factor is between 0 and 1
    R_noisy = clamp.(R_noisy, 0.0, 1.0)

    # Create a DataFrame
    df = DataFrame(Time=t_exprmnt, Recovery=R_noisy .* 100)  # Recovery as percentage

    # Save to CSV
    filename = filename * ".csv"
    CSV.write(filename, df)
end

end