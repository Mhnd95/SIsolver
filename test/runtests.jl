using SIsolver
using Test
using CSV
using DataFrames

@testset "SIsolver.jl" begin
    SIsolver.hello()
    
    # Test function for plot_error_vs_time_across_files with plot saving
    function test_plot_error_vs_time_across_files()
        # Create synthetic data for testing
        filenames = ["test_data1.csv", "test_data2.csv"]
        
        # Write synthetic data to CSV files
        data1 = DataFrame(Time = 1:10, Recovery = 0.1:0.1:1.0)
        data2 = DataFrame(Time = 1:10, Recovery = 0.2:0.1:1.1)
        CSV.write(filenames[1], data1)
        CSV.write(filenames[2], data2)
        
        initial_params = [0.33, 0.33, 0.33, 0.01, 0.01, 0.01, 0.0]
        save_path = "test_error_vs_time_across_files.png"
        
        # Call the function to test
        try
            SIsolver.plot_error_vs_time_across_files(filenames, initial_params, save_path)
            
            # Check if the file was created
            @test isfile(save_path) == true
            
            println("Test passed: Plot saved successfully.")
        catch e
            println("Test failed with error: $e")
            @test false
        end
    end
    
    # Run the test function
    test_plot_error_vs_time_across_files()
end