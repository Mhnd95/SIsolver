# SIsolver

Project to apply an algorithm that takes spontaneous imbibition data and predicts contact angle or any other missing parameter.  

Using the following sequence in Visual Studios:

file path:

Visual Studio = "data\\x.csv"
Github Codespace = "data/x.csv"

Initial optimization parameters:

θ = 0 radians
a[1:3] = 0.33
λ[1:3] = 0.01

julia --project=.

using SIsolver

SIsolver.optimize_parameters("data/0 hr.csv",[0.33,0.33,0.33,0.01,0.01,0.01,0.0])

SIsolver.optimize_params_across_files(["data/0 hr.csv","data/1 hr.csv"], [0.33,0.33,0.33,0.01,0.01,0.01,0.0])

Current step:

Status: The code takes one file and predicts all parameters from that input.
Goal: Need to take all data files and have the same "a" and "lambda" parameters across the different tests but independent "theta" for each case.

Trial:

SIsolver.plot_results(" hr.csv","trial")
