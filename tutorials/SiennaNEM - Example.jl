using SiennaNEM

using PowerSimulations

using Dates
using HiGHS

# input variables parameters
schedule_name = "schedule-1w"
scenario_name = 2

# data and system
data = SiennaNEM.get_data(
    "data/arrow",
    joinpath("data/arrow", schedule_name),
)
sys = SiennaNEM.create_system!(data)
SiennaNEM.add_ts!(
    sys, data;
    horizon=Hour(24),  # horizon of each time slice
    interval=Hour(24),  # interval between each time slice step in rolling horizon
    scenario_name=scenario_name,  # scenario number
)

# simulation
template_uc = SiennaNEM.build_problem_base_uc()
results = SiennaNEM.run_decision_model_loop(
    template_uc, sys;
    simulation_folder="examples/result/simulation_folder",
    simulation_name="$(schedule_name)_scenario-$(scenario_name)",
    decision_model_kwargs=(
        optimizer=optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01),
    ),
)