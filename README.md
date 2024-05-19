# SIsolver

Project to apply an algorithm that takes spontaneous imbibition data and predicts contact angle or any other missing parameter.  

SIsolver is a Julia module for optimizing parameters and visualizing data from spontaneous imbibition experiments. It predicts contact angles or other missing parameters using provided data.

## Table of Contents

- [Description](#description)
- [Installation](#installation)
- [Usage](#usage)
- [Current Status](#current-status)
- [Trial](#trial)
- [Results](#results)
  - [LBFGS Results](#lbfgs-results)
  - [NelderMead Results](#neldermead-results)
- [License](#license)
- [Contact](#contact)

## Description

SIsolver applies an algorithm to spontaneous imbibition data to predict contact angles or other missing parameters. It includes functions for reading data files, calculating shape factors, optimizing parameters, and generating various plots.

## Installation

To install SIsolver, clone the repository and add it to your Julia environment:

```sh
git clone https://github.com/Mhnd95/SIsolver.git
cd SIsolver
julia --project=.
```

In the Julia REPL, activate the environment and instantiate the dependencies:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Usage

Here's a basic example of how to use SIsolver:

```julia
using SIsolver

# Define the folder containing the data files and the file pattern
file_pattern = " hr.csv"

# Define the base save path for the plots
save_path = "results"

# Call the function to generate and save the plots and display the table
SIsolver.plot_results(file_pattern, save_path, max_iter=1000)
```

To optimize parameters for a single file:

```julia
using SIsolver

filename = "data/0 hr.csv"
optimized_params = SIsolver.optimize_parameters(filename)
```

## Current Status

* **Status**: The code runs smoothly and produces favorable results.
* **Goal**: Increase data set size and generate synthetic sets to test for code validity.

## Trial

```julia
using SIsolver

file_pattern = " hr.csv"
save_path = "trial"

SIsolver.plot_results(file_pattern, save_path, max_iter=1000)
```

## Results

### LBFGS Results

Optimized Parameters:

|File|a1|a2|a3|λ1|λ2|λ3|θ_deg
|:--------------:|:-----------|:-----------|:----------|:------------|:------------|:-----------|:---------
|data/0 hr.csv  | 0.00327045 | 0.00327045 |  0.967117 |  2.2942e-18 |  2.2942e-18 |  0.0129114 |  85.9344
|data/1 hr.csv  | 0.00327045 | 0.00327045 |  0.967117 |  2.2942e-18 |  2.2942e-18 |  0.0129114 |  87.5097
|data/6 hr.csv  | 0.00327045 | 0.00327045 |  0.967117 |  2.2942e-18 |  2.2942e-18 |  0.0129114 |  88.385
|data/12 hr.csv | 0.00327045 | 0.00327045 |  0.967117 |  2.2942e-18 |  2.2942e-18 |  0.0129114 |  88.7134
|data/24 hr.csv | 0.00327045 | 0.00327045 |  0.967117 |  2.2942e-18 |  2.2942e-18 |  0.0129114 |  89.0447
|data/48 hr.csv | 0.00327045 | 0.00327045 |  0.967117 |  2.2942e-18 |  2.2942e-18 |  0.0129114 |  89.2028

### NelderMead Results

Optimized Parameters:

|File|a1|a2|a3|λ1|λ2|λ3|θ_deg
|:--------------:|:-----------|:--------|:----------|:------------|:-----------|:-------------|:---------
|data/0 hr.csv  | 0.366921   | 7797.25 |  0.604348 |  0.00089564 |  1.16332e5 |  0.000895679 |  6.22454
|data/1 hr.csv  | 0.366921   | 7797.25 |  0.604348 |  0.00089564 |  1.16332e5 |  0.000895679 |  52.131
|data/6 hr.csv  | 0.366921   | 7797.25 |  0.604348 |  0.00089564 |  1.16332e5 |  0.000895679 |  66.5012
|data/12 hr.csv | 0.366921   | 7797.25 |  0.604348 |  0.00089564 |  1.16332e5 |  0.000895679 |  71.4522
|data/24 hr.csv | 0.366921   | 7797.25 |  0.604348 |  0.00089564 |  1.16332e5 |  0.000895679 |  76.3207
|data/48 hr.csv | 0.366921   | 7797.25 |  0.604348 |  0.00089564 |  1.16332e5 |  0.000895679 |  78.595

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

For questions or feedback, please contact me directly.

--------------------

### Terminal copied results

```sh
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
```
