# SIsolver
Project to apply an algorithm that takes spontaneous imbibition data and predicts contact angle or any other missing parameter.  

Command line running sequence:

SIsolver.read_data_file("data\\data.csv")

SIsolver.shape_factor(0.0,[0.0,0.0],[0.0,0.0])

SIsolver.experimental_results("data\\data.csv", 0.0, 0.238, 0.00, 0.0,0.0,0.0) #where 0.238 is the output of the previous line.


[![Build Status](https://github.com/Muhannad Alabdullateef/SIsolver.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/Muhannad Alabdullateef/SIsolver.jl/actions/workflows/CI.yml?query=branch%3Amaster)
