abstract type MatpowerRTEROPFSimpleInput <: AbstractInput end
using MAT

"""
    read_input(input_type::T, instance_path::String) where T<:Type{MatpowerInput}

Read instance in `instance_path` depending on `input_type`.\n
Return a structure OPFProblems.
"""
function read_input(input_type::T, instance_path::String) where T<:Type{MatpowerRTEROPFSimpleInput}
  file = matopen(instance_path)
  data = read(file, "mpc")
  close(file)
  # DataStructure and Gridstructure data:
  bus = SortedDict{String, SortedDict{String, Any}}()
  link = SortedDict{Link, SortedDict{String, Any}}()

  bus_id_line=SortedDict{Int, Int}()
  bus_id_name=SortedDict{Int, String}()

  baseMVA = data["baseMVA"]

  ## Building bus load and shunt information
  data_bus = data["bus"]
  bustype = data_bus[1:size(data_bus,1), 2]       # weather bus is 1:"PQ" (generators are not accounted for) or 2:"PV"
  for i in 1:size(data_bus,1)
      busname = bus_name(i)
      ## Adding MatpowerVolt structure (for each bus)
      if !haskey(bus, busname)
          bus[busname] = SortedDict{String, Any}()
      end
      bus[busname][volt_name()] = MatpowerVolt(busname, i, data_bus[i,13], data_bus[i,12])

      ## Matpower Load
      load = data_bus[i,3] + im*data_bus[i,4]
      if load != 0
        bus[busname]["Load"] = MatpowerLoad(load)
      end

      ## Matpower Shunt
      shunt = data_bus[i,5] + im*data_bus[i,6]
      if shunt != 0
        bus[busname]["Shunt"] = MatpowerROPFShunt(shunt)
      end

      bus_id_line[data_bus[i,1]] = i
      bus_id_name[data_bus[i,1]] = busname
  end

  ## Adding bus generator information
  gen2bus = SortedDict{Int, Int}()
  line2busgen = SortedDict{Int, Tuple{Int, Int}}()

  network, instance = split(instance_path, "\\")[end-1:end]
  gen_data =  data["gen"]
  instance = split(instance_path, "\\")[end][1:end-2]
  dict_Pinput = Dict{String, Tuple{Float64, Float64, Float64}}()
    genind, prevgen = 0, 0
  for i=1:size(gen_data,1)
    gen2bus[i] = bus_id_line[gen_data[i,1]]
    busname = bus_id_name[gen_data[i,1]]

    Pinput = 0
    S_min = S_max = 0
    if gen_data[i, 1] == prevgen
      genind += 1
      S_min = bus[busname]["Gen"].power_min
      S_max = bus[busname]["Gen"].power_max
      Pinput = bus[busname]["Gen"].active_power_input
    else
      prevgen = gen_data[i, 1]
      genind = 1
    end
    line2busgen[i] = (gen_data[i, 1], genind)

    if gen_data[i, 8] > 0 #generator is on
      S_min += gen_data[i,10] + im*gen_data[i,5] #Smin = Pmin + i Qmin
      S_max += gen_data[i,9] + im*gen_data[i,4] #Smax = Pmax + i Qmax
      Pinput += gen_data[i,2]
    end
    if bustype[bus_id_line[gen_data[i,1]]] == 3.0
      # println("reference bus : $i")
      network, instance = split(instance_path, "\\")[end-1:end]
      f = open("D:\\repo\\Test_Balthazar\\data_ROPF\\$(network)-$(instance[1:end-4])_refbus.csv", "w")
      write(f, "Gen_Sgen_$(bus_id_line[gen_data[i, 1]])_Re\n")
      close(f)
    end
    dict_Pinput["Gen_Sgen_$(bus_id_line[gen_data[i, 1]])_Re"] = (Pinput, real(S_min), real(S_max))
    bus[busname]["Gen"] = RTEMatpowerGenerator(S_min, S_max, Pinput, true)
  end
  network, instance = split(instance_path, "\\")[end-1:end]
  f = open("D:\\repo\\Test_Balthazar\\data_ROPF\\$(network)-$(instance[1:end-4])_generation.csv", "w")
  for (gen_name, tuple) in dict_Pinput
    Pinput = tuple[1]
    Pmin = tuple[2]
    Pmax = tuple[3]
    write(f, "$(gen_name) ; $Pinput ; ; $(Pmin) ; $(Pmax) \n")
  end
  close(f)

  ## building link information
  data_branch = data["branch"]
  for i=1:size(data_branch,1)
    linkname = Link(bus_id_name[data_branch[i,1]], bus_id_name[data_branch[i,2]])
    if data_branch[i, 11] == 0
      println("$link $(linkname.orig)⟶$(linkname.dest) breaker out of service !")
    else
      rs, xs, bc = data_branch[i,3:5]
      Smax = data_branch[i,6]
      τ, θ = data_branch[i, 9:10]
      if !haskey(link, linkname)
        link[linkname] = SortedDict{String, Any}()
      end

      nb_elem = length(link[linkname])
      link[linkname]["LinkMP_$(nb_elem+1)"] = MatpowerLine_π(baseMVA, rs, xs, bc, τ, θ)
    end
  end


  bus_line_id = SortedDict([val=>key for (key,val) in bus_id_line])

  ## Removing generators from PQ buses
  for (busname, bus_elems) in bus
    i = bus_elems["Volt"].busid
    if bustype[i] == 1
      for (elemid, elem) in bus_elems
        if typeof(elem) == MatpowerGenerator
          delete!(bus_elems, elemid)
        end
      end
    end
  end

  ## Adding null generator to reference bus
  refind = findall(bustype .== 3)
  length(refind) == 1 || warn("/! refind = $refind")
  bus["BUS_$(refind[1])"]["Gen_reference"] = MatpowerGenerator(0,0,3,[0 0 0],true)

  ds = DataSource(bus, link)

  ## Building GridStructure
  node_linksin = node_linksout = SortedDict{String, SortedSet{Link}}()
  node_vars = SortedDict{String, SortedDict{String, Variable}}()
  link_vars = SortedDict{Link, SortedDict{String, Variable}}()
  gs = GridStructure(basecase_scenario_name(), node_linksin, node_linksout)

  node_formulations = SortedDict{String, SortedDict{String, Symbol}}()
  link_formulations = SortedDict{Link, SortedDict{String, Symbol}}()
  mp = MathematicalProgramming(node_formulations, link_formulations, node_vars, link_vars)
  return OPFProblems(basecase_scenario_name()=>Scenario(ds, gs, mp))
end
