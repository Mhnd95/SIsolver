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

Status: The code runs smoothly and produces favorable results.
Goal: Increase data set size and generate synthetic sets to test for code validity.

Trial:

SIsolver.plot_results(" hr.csv","trial")

LBFG results:

Optimized Parameters:
6×8 DataFrame
 Row │ File            a1          a2          a3        λ1          λ2          λ3         θ_deg
     │ String          Float64     Float64     Float64   Float64     Float64     Float64    Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────────
   1 │ data/0 hr.csv   0.00327045  0.00327045  0.967117  2.2942e-18  2.2942e-18  0.0129114  85.9344
   2 │ data/1 hr.csv   0.00327045  0.00327045  0.967117  2.2942e-18  2.2942e-18  0.0129114  87.5097
   3 │ data/12 hr.csv  0.00327045  0.00327045  0.967117  2.2942e-18  2.2942e-18  0.0129114  88.7134
   4 │ data/24 hr.csv  0.00327045  0.00327045  0.967117  2.2942e-18  2.2942e-18  0.0129114  89.0447
   5 │ data/48 hr.csv  0.00327045  0.00327045  0.967117  2.2942e-18  2.2942e-18  0.0129114  89.2028
   6 │ data/6 hr.csv   0.00327045  0.00327045  0.967117  2.2942e-18  2.2942e-18  0.0129114  88.385

NelderMead results:

Optimized Parameters:
6×8 DataFrame
 Row │ File            a1        a2       a3        λ1          λ2         λ3           θ_deg
     │ String          Float64   Float64  Float64   Float64     Float64    Float64      Float64  
─────┼───────────────────────────────────────────────────────────────────────────────────────────
   1 │ data/0 hr.csv   0.366921  7797.25  0.604348  0.00089564  1.16332e5  0.000895679   6.22454
   2 │ data/1 hr.csv   0.366921  7797.25  0.604348  0.00089564  1.16332e5  0.000895679  52.131
   3 │ data/12 hr.csv  0.366921  7797.25  0.604348  0.00089564  1.16332e5  0.000895679  71.4522
   4 │ data/24 hr.csv  0.366921  7797.25  0.604348  0.00089564  1.16332e5  0.000895679  76.3207
   5 │ data/48 hr.csv  0.366921  7797.25  0.604348  0.00089564  1.16332e5  0.000895679  78.595
   6 │ data/6 hr.csv   0.366921  7797.25  0.604348  0.00089564  1.16332e5  0.000895679  66.5012
