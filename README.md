# SIsolver

Project to apply an algorithm that takes spontaneous imbibition data and predicts contact angle or any other missing parameter.  

Using the following sequence in Visual Studios:

file path = "data\\data.csv"

Initial optimization parameters:

θ = 0 radians
a[1:3] = 0.33
λ[1:3] = 0.01

julia --project=.

using SIsolver

SIsolver.optimize_parameters("data\\data.csv",[0.33,0.33,0.33,0.01,0.01,0.01,0.0])

Current step:

* The code takes one file and predicts all parameters from that input. Need to take all data files and have the same a and \lambda across the different tests but independent \theta for each case.

Trials:

SIsolver.optimize_params_across_files(["data\\data.csv","data\\1 hr.csv"], [0.33,0.33,0.33,0.01,0.01,0.01,0.0])

[![Build Status](https://github.com/Muhannad Alabdullateef/SIsolver.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Muhannad Alabdullateef/SIsolver.jl/actions/workflows/CI.yml?query=branch%3Amain)
