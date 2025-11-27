

#%% Initialise all dependencies and packages
using Pkg
Pkg.develop(path="../PRASNEM.jl") # Adjust path as needed
Pkg.develop(path="../PISP.jl") # Adjust path as needed
Pkg.develop(path="../SchedNEM.jl") # Adjust path as needed
using SchedNEM
using PISP
using PRASNEM
using PRAS
using Plots
using CSV
using DataFrames
using Dates
using Gurobi

#%% ================= Modelling storage operation details ===============================

input_folder = joinpath(pwd(), "..", "data", "nem12")
output_folder = joinpath(pwd(), "..", "data", "nem12", "pras")
tyear = 2030

start_dt = DateTime("$(tyear)-01-01 00:00:00", dateformat"yyyy-mm-dd HH:MM:SS")
end_dt = DateTime("$(tyear)-01-14 23:00:00", dateformat"yyyy-mm-dd HH:MM:SS")
timeseries_folder = joinpath(input_folder, "schedule-$(tyear)")
sys = PRASNEM.create_pras_system(start_dt, end_dt, input_folder, timeseries_folder; output_folder=output_folder, gentech_excluded=["DSP"])
sys.regions.load[:,:] .= round.(Int, sys.regions.load[:,:] * 1.25)
sys.attrs["case"] = sys.attrs["case"] * "_load+25pct"

m = SchedNEM.build_operation_model(sys; optimisation_window=48, move_forward=24, input_folder=input_folder, optimiser=Gurobi.Optimizer())
res = SchedNEM.run_operation_model(m, sys)

#%% ========================================

# AEMO energy capacity derating
sys_stor_derated = deepcopy(sys)
sys_stor_derated = SchedNEM.updateEnergyDerating(sys_stor_derated) # Res not needed, as derating is fixed

# Compare with expecation dispatch
sys_stor_fixed = deepcopy(sys)
sys_stor_fixed = SchedNEM.updateMarketExpectationDispatch(sys_stor_fixed, res; include_genstorage=false)

#%% Run PRAS with different storage modelling approaches
simspecs = SequentialMonteCarlo(samples=200, seed=110);
resultspecs = (Shortfall(), StorageEnergy(),);

sf_greedy, se_greedy, = assess(sys, simspecs, resultspecs...)
sf_derated, se_derated,  = assess(sys_stor_derated, simspecs, resultspecs...)
sf_expectation, se_expectation,  = assess(sys_stor_fixed, simspecs, resultspecs...)

#%% ========================================
# Display results


plot(sum(se_greedy.energy_mean, dims=1)[:] ./ 66, label="Greedy", legend=:topright, 
    xlabel="Time [hours]", ylabel="Aggregated State of Charge", linewidth=2, dpi=300,
    size=(400,400), yticks=false)
plot!(sum(se_derated.energy_mean, dims=1)[:] ./ 66, label="Energy derated", linewidth=2)
plot!(sum(res.stor_energy, dims=1)[:] ./ 66, label="Market-based dispatch", linewidth=2)
ylims!(0, 550)
xlims!(0, 200)
savefig("./scripts/AR-PST Monthly Meetings/figures/storage_operation_comparison.png")

