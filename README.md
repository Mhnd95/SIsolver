# SIsolver
Project to apply an algorithm that takes spontaneous imbibition data and predicts contact angle or any other missing parameter.  

Using the following sequence in Visual Studios:

file path = "data\\data.csv"

julia --project=.

using SIsolver

SIsolver.R_calculated("data\\data.csv", 0.0, [0.33,0.33,0.33], [0.01,0.01,0.01])

Current issue:

* fix the R_calcualted as exp function cannot take vector as input, the vector in this case is t_experimental

[![Build Status](https://github.com/Muhannad Alabdullateef/SIsolver.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Muhannad Alabdullateef/SIsolver.jl/actions/workflows/CI.yml?query=branch%3Amaster)
