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

SIsolver.global_objective("data\\data.csv", 0.0, [0.33,0.33,0.33], [0.01,0.01,0.01])


Current step:

* Work on optimization code functions

[![Build Status](https://github.com/Muhannad Alabdullateef/SIsolver.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Muhannad Alabdullateef/SIsolver.jl/actions/workflows/CI.yml?query=branch%3Amain)
