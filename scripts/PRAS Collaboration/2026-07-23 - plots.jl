


using PRASNEM
using PRAS
using Dates
using DataFrames
using Statistics
using CSV
using Plots

include(joinpath(@__DIR__, "utils.jl"));

#%%

# Create system for 2035 with reference year 2011 and 2019, with high and low demand
sys11 = create_system(2035, 2011, 10)
sys19 = create_system(2035, 2019, 50)
show_system_details(sys11)
show_system_details(sys19)


#%%

# Run adequacy assessment with 1000 samples, with seed 123

run_base_adequacy_assessment(sys11)
run_base_adequacy_assessment(sys19)

#%%

# Collect event details for 1000 samples, with seed 123 and following

ev_details_11 = get_event_details(sys11)
ev_details_19 = get_event_details(sys19)

ev_details_11.annual_demand_mwh .= sum(sys11.regions.load)
ev_details_19.annual_demand_mwh .= sum(sys19.regions.load)

e = vcat(ev_details_11, ev_details_19)

CSV.write(joinpath(@__DIR__, "event_details_2035_2011.csv"), ev_details_11)

e = ev_details_11

#%% Risk profile plots

xlims_overall = (0, 18) #(0, maximum(r.events.duration_hrs) + 1)
ylims_overall = (-2, 150) # maximum(r.events.magnitude_mwh ./ r.events.total_load_mw) .* 1e3 + 0.001)

yticks_overall = 0:30:ylims_overall[2]
xticks_overall = 0:3:xlims_overall[2]

kwargs_scatter = (label="", markersize=3, xlabel= "Duration [hrs]",
    xlims=xlims_overall, ylims=ylims_overall, ylabel="Event ENS\n[ppm of annual demand]",
   yticks=yticks_overall, xticks=xticks_overall, legend=false,
   size=(350,300), dpi=500)

scatter(e.length, e.sum ./ e.annual_demand_mwh .* 1e6, c=1; kwargs_scatter...)
#savefig(joinpath(@__DIR__, "figures","risk_profile_2035_2011.png"))

stephist(e.length, bins=0:1:xlims_overall[2], label="", 
   lw=2, c=1, xlabel="Duration [hrs]", ylabel="Number of events [#/yr]",
   xlims=xlims_overall, xticks=xticks_overall, ylims=(0, 12000),
   yticks=(0:1000:12000, string.(collect(0:1000:12000) ./ 10000)),
   size=(600,300), dpi=500)
savefig(joinpath(@__DIR__, "figures","stephist_duration.png"))

stephist(e.sum ./ e.annual_demand_mwh .* 1e6, bins=0:5:ylims_overall[2], label="", 
   lw=2, c=1, xlabel="Event ENS [ppm of annual demand]", ylabel="Number of events [#/yr]",
   #xlims=xlims_overall, xticks=xticks_overall, ylims=(0, 12000),
   xticks=(0:30:150), xlims=(-2, 150), ylims=(0, 25000),
   yticks=(0:5000:50000, string.(collect(0:5000:50000) ./ 10000)),
   size=(600,300), dpi=500)
savefig(joinpath(@__DIR__, "figures","stephist_ens.png"))





#%% Time of day plot

bins = 0:1:24

critical_hours_start = zeros(24)

for row in eachrow(e)
   critical_hours_start[round(Int, row.start_index) % 24 + 1] += 1
end
critical_hours_start ./= nrow(e)

x = repeat(0:1:23, inner=2)
y = vcat(zeros(1), repeat(critical_hours_start, inner=2)[1:end-1])
plot(x, y * 100, label="", lw=2, c= 1, fillalpha=0.5, fillcolor=[1],
   yticks=(0:20:100), ylims=(-1, 100), ylabel="Probability of event start [%]",
   size=(400,300), dpi=500, xlabel="Hour of day [hr]", xticks=0:3:24)
savefig(joinpath(@__DIR__, "figures","event_start_time.png"))

#%% Month of year plot 

bins = 0:1:12
critical_month_start = zeros(12)

for row in eachrow(e)
   critical_month_start[month(start_dt + Hour(row.start_index))] += 1
end
critical_month_start ./= nrow(e)

x = repeat(1:1:12, inner=2)
y = vcat(zeros(1), repeat(critical_month_start, inner=2)[1:end-1])
plot(x, y * 100, label="", lw=2, c= 1, fillalpha=0.5, fillcolor=[1],
   yticks=(0:20:100), ylims=(-1, 100), ylabel="Probability of event start [%]",
   size=(400,300), dpi=500, xlabel="Month of year", xticks=(collect(1:1:12) .+ 0.5, Dates.monthabbr.(1:12)))
savefig(joinpath(@__DIR__, "figures","event_start_month.png"))

