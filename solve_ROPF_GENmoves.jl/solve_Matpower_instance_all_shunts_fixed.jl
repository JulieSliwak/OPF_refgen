include("generic_functions.jl")
#example
output_instance_path = "..\\data_ROPF"
output_decomposition_path = "..\\data_sdp"
list_generation = ["generation"]
generation_files_path = "..\\data_ROPF"
network_name = "118nodes"

networks = ["118nodes_disconnect"]

for network_name in networks
    f = open("results_test_allshunts1_$(network_name).csv", "w")
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
            UB = Inf
            UB, LB = solve1_refgen_allshunts1(ROPF)
            println("UB=$UB ; LB=$LB")
            if UB != Inf
                cp("solutions\\solve1_allshunts1_refgen_solution_$(instance_name).csv", "solutions_mp_data\\allshunts1_refgen_solution_$(instance_name).csv")
            end
            gap = (UB-LB)/UB
            write(f, "$numero ; $UB ; $LB ; $gap \n")
        end
    end
    close(f)
end
