include("generic_functions.jl")
#example
output_instance_path = "..\\data_ROPF"
output_decomposition_path = "..\\data_sdp"
list_generation = ["generation"]
generation_files_path = "..\\data_ROPF"
max_time = 600 #10 minutes
network_name = "14nodes"

f = open("results_test_$(network_name)_v+-0.02_BB$(max_time).csv", "w")
write(f, "Numero ; UB ; LB ; Gap (%)\n")
for i in 0:2499
    # network_name, numero = split(instance_name, '-')
    numero = "test_$i"
    instance_name = network_name*"-"*numero
    matpower_instance_path = "..\\mp_data\\$(network_name)\\$(numero).mat"
    for generation in list_generation
        construct_dat_file_ROPF(instance_name, matpower_instance_path, output_instance_path)
        # generate_clique_decomposition(instance_name, matpower_instance_path, output_decomposition_path)
        ROPF = ROPF_infos(instance_name,
        matpower_instance_path,
        output_instance_path,
        "cholesky",
        output_decomposition_path,
        generation,
        generation_files_path)
        UB, LB = solve1_refgen(ROPF)
        println("UB=$UB ; LB=$LB")
        UB_BB, LB = solve2_refgen(ROPF, max_time)
        println("After B&B")
        println("UB=$UB_BB ; LB=$LB")
        if UB_BB != Inf
            if UB_BB < UB
                cp("solutions\\BB_refgen_solution_$(instance_name).csv", "solutions_mp_data\\refgen_solution_$(instance_name).csv")
            else
                cp("solutions\\solve1_refgen_solution_$(instance_name).csv", "solutions_mp_data\\refgen_solution_$(instance_name).csv")
            end
        end
        gap = (min(UB, UB_BB)-LB)/(min(UB, UB_BB))
        write(f, "$numero ; $(min(UB, UB_BB)) ; $LB ; $gap \n")
    end
end
close(f)
