include(joinpath("..","src", "PowSysMod_body.jl"))

function build_matpower_instance(instance_path, typeofinput)
    #load and build complex problem
    OPFpbs = load_OPFproblems(typeofinput, instance_path)
    pb_global = build_globalpb!(OPFpbs)
    #convert complex problem to real problem
    pb_global_real = pb_cplx2real(pb_global)
end

#typeofinput = MatpowerSimpleInput
#typeofinput = MatpowerROPFSimpleImaxInput

#problem can be exported using function 'export_to_dat'
# export_to_dat(pb_global_real, output_path, filename="$instance.dat")
