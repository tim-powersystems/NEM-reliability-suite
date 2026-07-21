

function create_system(target_year, ref_trace, poe; DERcase="base", scenario=2, base_path="Z:/")
    add_lines = PRASNEM.get_added_lines_per_year()
    DER_parameters = PRASNEM.get_DER_parameters(; case=DERcase)

    start_dt = DateTime("$(target_year)-01-01 00:00:00", dateformat"yyyy-mm-dd HH:MM:SS")
    end_dt = DateTime("$(target_year)-12-31 23:00:00", dateformat"yyyy-mm-dd HH:MM:SS")

    base_folder_pisp = joinpath(base_path, "pisp-datasets", "base")
    base_folder_pras = joinpath(base_path, "pras-files", DERcase)
    pisp_input_folder = joinpath(base_folder_pisp, "out-ref$(ref_trace)-poe$(poe)", "csv")
    timeseries_folder = "schedule-$(target_year)"
    pras_folder = joinpath(base_folder_pras, "out-ref$(ref_trace)-poe$(poe)",)

    sys = PRASNEM.create_pras_system(start_dt, end_dt, pisp_input_folder, timeseries_folder; 
                               line_alias_included=add_lines[target_year], output_folder=pras_folder,
                               DER_parameters=DER_parameters, scenario=scenario)
    
    return sys
end


function show_system_details(sys)
   # Print system details
   println("System details:")
   println("Peak demand: ", maximum(sum(sys.regions.load, dims=1)), " MW")

   idxs_vre = findall(x -> x in ["RoofPV", "LargePV", "Wind"], sys.generators.categories)
   idxs_thermal = findall(x -> !(x in ["RoofPV", "LargePV", "Wind"]), sys.generators.categories)
   idxs_hydro = findall(x -> !(x in ["PS"]), sys.generatorstorages.categories)
   idxs_phes = findall(x -> x in ["PS"], sys.generatorstorages.categories)
   idxs_vre_per_region = [intersect(idxs_vre, sys.region_gen_idxs[region]) for region in 1:12]
   idxs_thermal_per_region = [intersect(idxs_thermal, sys.region_gen_idxs[region]) for region in 1:12]

   println("Total VRE capacity: ", sum(maximum(sys.generators.capacity[idxs_vre,:], dims=2)), " MW")
   println("Total thermal capacity: ", sum(maximum(sys.generators.capacity[idxs_thermal,:], dims=2)), " MW")
   println("Total hydro capacity: ", sum(maximum(sys.generatorstorages.discharge_capacity[idxs_hydro,:], dims=2)), " MW")
   println("Total storage capacity: ", sum(maximum(sys.generatorstorages.discharge_capacity[idxs_phes,:], dims=2)) + sum(maximum(sys.storages.discharge_capacity, dims=2)), " MW")

   annual_cf_per_region = [sum(sys.generators.capacity[idxs_vre_per_region[region],:]) / sum(maximum(sys.generators.capacity[idxs_vre_per_region[region],:], dims=2)) for region in 1:12] ./ 8760 .* 100
   println("Annual capacity factor: ", round(mean(annual_cf_per_region), digits=2), " % (min: ", round(minimum(annual_cf_per_region), digits=2), ", max: ", round(maximum(annual_cf_per_region), digits=2), ")")

end

function run_base_adequacy_assessment(sys; samples=1000, seed=123)
   println("Running adequacy assessment with $(samples) samples, seed=$(seed)")
   sf, = assess(sys, SequentialMonteCarlo(samples=samples, seed=seed), Shortfall())
   println("LOLE: ", LOLE(sf))
   println("EUE: ", EUE(sf))
   println("NEUE: ", NEUE(sf))
end

function get_event_details(sys; samples=1000, seed=123)
   all_events = DataFrame()

   n_runs = round(Int, samples / 100)

   Threads.@threads for i in 1:n_runs
      println("$i / $n_runs: Collecting event details for 100 samples, seed=$(seed)")
      sfsamples, = assess(sys, SequentialMonteCarlo(samples=samples, seed=seed), ShortfallSamples())
      df_events = PRASNEM.get_all_event_details(sfsamples)
      append!(all_events, df_events)
      seed = seed + 1
   end

   return all_events
end