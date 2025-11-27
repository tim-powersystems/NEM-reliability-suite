
#%% Initialise all dependencies and packages
using Pkg
Pkg.develop(path="../PRASNEM.jl") # Adjust path as needed
Pkg.develop(path="../PISP.jl") # Adjust path as needed
using PISP
using PRASNEM
using PRAS
using Plots
using CSV
using DataFrames
using Dates

#%% ================= Modelling demand flexibility ===============================

simspec = SequentialMonteCarlo(samples=200, seed=110);
respec = (Shortfall(), DemandResponseEnergy(),);
tyear = 2030

sf_res_dsp = zeros(Float64, 2, 1) # Store NEUE results for 4 different storage dispatches, 11 target years

added_lines_per_year = PRASNEM.get_added_lines_per_year(2) # Scenario 2: ISP Step Change

input_folder = joinpath(pwd(), "..", "data", "nem12")
output_folder = joinpath(pwd(), "..", "data", "nem12", "pras")


# Get the original system first
println("=============================== $tyear ==============================")
start_dt = DateTime("$(tyear)-01-01 00:00:00", dateformat"yyyy-mm-dd HH:MM:SS")
end_dt = DateTime("$(tyear)-12-31 23:00:00", dateformat"yyyy-mm-dd HH:MM:SS")
timeseries_folder = joinpath(input_folder, "schedule-$(tyear)")
sys = PRASNEM.create_pras_system(start_dt, end_dt, input_folder, timeseries_folder; output_folder=output_folder)
sys_without_dsp = PRASNEM.create_pras_system(start_dt, end_dt, input_folder, timeseries_folder; output_folder=output_folder, gentech_excluded=["DSP"])

#%% Run PRAS with and without DSP
sf_with_dsp, dr_energy = assess(sys, simspec, respec...);
sf_without_dsp, _ = assess(sys_without_dsp, simspec, respec...);

println("Year $tyear: Shortfall with DSP = $(NEUE(sf_with_dsp)) MWh, without DSP = $(NEUE(sf_without_dsp)) MWh")
println("Without DSP:")
PRASNEM.NEUE_area(sys_without_dsp, sf_without_dsp; bus_file_path=joinpath(input_folder, "Bus.csv"));
println("With DSP:")
PRASNEM.NEUE_area(sys, sf_with_dsp; bus_file_path=joinpath(input_folder, "Bus.csv"));
