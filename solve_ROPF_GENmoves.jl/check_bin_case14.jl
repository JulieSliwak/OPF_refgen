function check_bin_1()
    nb_var_1 = 0
    nb_files = 0
    for file in readdir("solutions_mp_data")
        nb_files +=1
        lines = readlines(joinpath("solutions_mp_data", file))
        for line in lines
            var, value = split(line, ';')
            if var == "Shunt_9"
                if parse(Float64, value) == 1.0
                    nb_var_1 +=1
                else
                    println(file)    
                end
            end
        end
    end
    return nb_files, nb_var_1
end

 nb_files, nb_var_1 = check_bin_1()
println("nb files = $nb_files")
println("nb var 1 = $nb_var_1")
