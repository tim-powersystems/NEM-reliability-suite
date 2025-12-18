using SiennaNEM
using PowerSimulations
using HiGHS
using Dates

system_data_dir = "data/arrow"
ts_data_dir     = joinpath(system_data_dir, "schedule-1w")
scenario        = 1
file_format     = "arrow"

data = read_system_data(system_data_dir; file_format=file_format);
read_ts_data!(data, ts_data_dir; file_format=file_format);
add_tsf_data!(data, scenario=scenario);
update_system_data_bound!(data);
clean_ts_data!(data)

sys = create_system!(data);
SiennaNEM.add_ts!(
    sys, data;
    horizon=Hour(24),
    interval=Hour(24),
    scenario=scenario,
)
template_uc = SiennaNEM.build_problem_base_uc();
solver = optimizer_with_attributes(HiGHS.Optimizer, "mip_rel_gap" => 0.01);

hours = Hour(24)
problem = DecisionModel(template_uc, sys; optimizer=solver, horizon=hours)
@time build!(problem; output_dir=mktempdir())
@time solve!(problem)
res = OptimizationProblemResults(problem)
