using MAT

file = matopen("D:\\repo\\Test_Balthazar\\mp_data\\118nodes\\test_0.mat")
varnames = names(file)
#println(varnames)
data = read(file, "mpc")
println(data["gen"][:,1:2])
# println(data["branch"][:,3:4])
close(file)
