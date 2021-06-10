abstract type MatpowerRTEROPFSimpleInput <: AbstractInput end


"""
    read_input(input_type::T, instance_path::String) where T<:Type{MatpowerInput}

Read instance in `instance_path` depending on `input_type`.\n
Return a structure OPFProblems.
"""
function read_input(input_type::T, instance_path::String) where T<:Type{MatpowerRTEROPFSimpleInput}
  data =  read(instance_path, "mpc")
  close(instance_path)
  # DataStructure and Gridstructure data:
  bus = SortedDict{String, SortedDict{String, Any}}()
  link = SortedDict{Link, SortedDict{String, Any}}()

  bus_id_line=SortedDict{Int, Int}()
  bus_id_name=SortedDict{Int, String}()

  checkfor(data, 2, "mpc.baseMVA")
  baseMVA = data[2,3]

  ## Building bus load and shunt information
  i_debut, i_fin = find_numarray(1, data)
  bustype = data[i_debut:i_fin, 2]                # weather bus is 1:"PQ" (generators are not accounted for) or 2:"PV"
  checkfor(data, i_debut-1, "mpc.bus")
  for i=i_debut:i_fin
    id = i-i_debut+1
    busname = bus_name(id)
    ## Adding MatpowerVolt structure (for each bus)
    if !haskey(bus, busname)
        bus[busname] = SortedDict{String, Any}()
    end
    bus[busname][volt_name()] = MatpowerVolt(busname, id, data[i,13], data[i,12])

    ## Matpower Load
    load = data[i,3] + im*data[i,4]
    if load != 0
      bus[busname]["Load"] = MatpowerLoad(load)
    end

    ## Matpower Shunt
    shunt = data[i,5] + im*data[i,6]
    if shunt != 0
      bus[busname]["Shunt"] = MatpowerROPFShunt(shunt)
    end

    bus_id_line[data[i,1]] = id
    bus_id_name[data[i,1]] = busname
  end

  ## Adding bus generator information
  gen2bus = SortedDict{Int, Int}()
  line2busgen = SortedDict{Int, Tuple{Int, Int}}()

  i_debut, i_fin = find_numarray(i_fin+1, data)
  checkfor(data, i_debut-1, "mpc.gen")
  instance = split(instance_path, "\\")[end][1:end-2]
  dict_Pinput = Dict{String, Tuple{Float64, Float64, Float64}}()
    genind, prevgen = 0, 0
  for i=i_debut:i_fin
    gen2bus[i-i_debut+1] = bus_id_line[data[i,1]]
    busname = bus_id_name[data[i,1]]

    Pinput = 0
    S_min = S_max = 0
    if data[i, 1] == prevgen
      genind += 1
      S_min = bus[busname]["Gen"].power_min
      S_max = bus[busname]["Gen"].power_max
      Pinput = bus[busname]["Gen"].active_power_input
    else
      prevgen = data[i, 1]
      genind = 1
    end
    line2busgen[i-i_debut+1] = (data[i, 1], genind)

    if data[i, 8] > 0 #generator is on
      S_min += data[i,10] + im*data[i,5] #Smin = Pmin + i Qmin
      S_max += data[i,9] + im*data[i,4] #Smax = Pmax + i Qmax
    end
    Pinput = 0.5*(real(S_min) + real(S_max))
    # Pinput += data[i,2]
    dict_Pinput["Gen_Sgen_$(bus_id_line[data[i, 1]])_Re"] = (Pinput, real(S_min), real(S_max))
    bus[busname]["Gen"] = RTEMatpowerGenerator(S_min, S_max, Pinput, true)
  end
  # f = open("D:\\repo\\Test_Balthazar\\data_ROPF\\$(instance)_generation1.csv", "w")
  # for (gen_name, tuple) in dict_Pinput
  #   Pinput = tuple[1]
  #   Pmin = tuple[2]
  #   Pmax = tuple[3]
  #   write(f, "$(gen_name) ; $Pinput ; ; $(Pmin) ; $(Pmax) \n")
  # end
  # close(f)

  ## building link information
  i_debut, i_fin = find_numarray(i_fin+1, data)
  checkfor(data, i_debut-1, "mpc.branch")
  for i=i_debut:i_fin
    linkname = Link(bus_id_name[data[i,1]], bus_id_name[data[i,2]])
    if data[i, 11] == 0
      #@warn("$link $(linkname.orig)⟶$(linkname.dest) breaker out of service !")
    else
      rs, xs, bc = data[i,3:5]
      Smax = data[i,6]
      τ, θ = data[i, 9:10]
      if !haskey(link, linkname)
        link[linkname] = SortedDict{String, Any}()
      end

      nb_elem = length(link[linkname])
      link[linkname]["LinkMP_$(nb_elem+1)"] = MatpowerLine_π(baseMVA, rs, xs, bc, τ, θ)
    end
  end

  i_debut, i_fin = find_numarray(i_fin+1, data)
  if data[i_debut-1, 1] == "mpc.areas"
    @warn("Not loading mpc.areas data.")
    i_debut = i_fin+1
    i_debut, i_fin = find_numarray(i_fin+1, data)
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


## Utils functions
function find_numarray(i_start, data)
  i_debut = i_start
  while !isa(data[i_debut, 1], Int)
    i_debut+=1
  end
  i_fin=i_debut
  while !isa(data[i_fin,1], SubString)
    i_fin += 1
  end
  i_debut, i_fin-1
end

checkfor(data, line_ind, name) = (data[line_ind, 1] == name) || error("Expected ", name, " at line ", line_ind, ", got ", data[line_ind,1], " instead.")

function load_matpower(filename)
  instance_name = split(filename, '.')[1]

  touch(instance_name*".temp")
  f = open(filename)
  out = open(instance_name*".temp", "w")

  # removing all ';' at end of lines
  while !eof(f)
    line = readline(f)
    if length(line) > 0 && line[1] != '%' && line[1] != 'f'
      s = split(line, ";")
      println(out, s[1])
    end
  end
  close(f)
  close(out)

  data = DelimitedFiles.readdlm(instance_name*".temp")
  rm(instance_name*".temp")
  return data
end
