using DelimitedFiles, DataStructures

function get_features(instance_path::String)
  data = load_matpower(instance_path)
  costs_generators = Dict{String, Tuple{Float64,Float64}}() # bus g => (cg, kg)
  loads = Dict{String, Complex}() # bus => Sl= Pl + im*Ql
  shunts = Dict{String, Complex}() # bus n => g - im *b
  bounds_voltage = Dict{String, Tuple{Float64, Float64}}() # bus => (vmin,vmax)
  bounds_realpower = Dict{String, Tuple{Float64, Float64}}() # bus => (Pmin,Pmax)
  bounds_imagpower = Dict{String, Tuple{Float64, Float64}}() # bus => (Qmin,Qmax)
  max_current = Dict{Tuple{String,String}, Float64}() # line=(bus,bus) => imax

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
    busname = "BUS_$id"
    bounds_voltage[busname] = (data[i,13], data[i,12])
    loads[busname] = data[i,3] + im*data[i,4]
    if data[i,5] == data[i,6] == 0
      #no shunts
    else
      shunts[busname] = data[i,5] - im*data[i,6]
    end
    bus_id_line[data[i,1]] = id
    bus_id_name[data[i,1]] = busname
  end

  ## Adding bus generator information
  gen2bus = SortedDict{Int, Int}()
  line2busgen = SortedDict{Int, Tuple{Int, Int}}()

  i_debut, i_fin = find_numarray(i_fin+1, data)
  checkfor(data, i_debut-1, "mpc.gen")

  genind, prevgen = 0, 0
  for i=i_debut:i_fin
    gen2bus[i-i_debut+1] = bus_id_line[data[i,1]]
    busname = bus_id_name[data[i,1]]

    S_min = S_max = 0
    if data[i, 1] == prevgen
      println("Several generators at bus $busname, $(data[i,1])")
      genind += 1
      (Pmin, Pmax) = bounds_realpower[busname]
      (Qmin, Qmax) = bounds_imagpower[busname]
      S_min = Pmin + im*Qmin
      S_max = Pmax + im*Qmax
      # S_min = bus[busname]["Gen"].power_min
      # S_max = bus[busname]["Gen"].power_max
    else
      prevgen = data[i, 1]
      genind = 1
    end
    line2busgen[i-i_debut+1] = (data[i, 1], genind)

    if data[i, 8] > 0 #generator is on
      S_min += data[i,10] + im*data[i,5] #Smin = Pmin + i Qmin
      S_max += data[i,9] + im*data[i,4] #Smax = Pmax + i Qmax
    end
    bounds_realpower[busname] = (real(S_min), real(S_max))
    bounds_imagpower[busname] = (imag(S_min), imag(S_max))

  end

  ## building link information
  i_debut, i_fin = find_numarray(i_fin+1, data)
  checkfor(data, i_debut-1, "mpc.branch")
  for i=i_debut:i_fin
    if data[i, 11] == 0
      #@warn("$link $(linkname.orig)⟶$(linkname.dest) breaker out of service !")
    else
      rs, xs, bc = data[i,3:5]
      Smax = data[i,6]
      τ, θ = data[i, 9:10]
      if Smax != 0
        max_current[(bus_id_name[data[i,1]], bus_id_name[data[i,2]])] = Smax/baseMVA
      end
    end
  end

  i_debut, i_fin = find_numarray(i_fin+1, data)
  if data[i_debut-1, 1] == "mpc.areas"
    @warn("Not loading mpc.areas data.")
    i_debut = i_fin+1
    i_debut, i_fin = find_numarray(i_fin+1, data)
  end

  bus_line_id = SortedDict([val=>key for (key,val) in bus_id_line])

  ## Adding generator cost information
  checkfor(data, i_debut-1, "mpc.gencost")
  genind, cur_genind = 0, data[i_debut, 1]
  for i=i_debut:i_fin
    buslineid, _ = line2busgen[i-i_debut+1]
    busid = bus_id_line[buslineid]

    busname = "BUS_$busid"

    cost_degree = data[i, 4]
    cost_coeffs = data[i, 5:(5+cost_degree-1)]
    if cost_degree == 2
      costs_generators[busname] = (cost_coeffs[1], cost_coeffs[2])
    elseif cost_degree == 3
      costs_generators[busname] = (cost_coeffs[2], cost_coeffs[3])
    end
  end
  return costs_generators, loads, shunts, bounds_voltage, bounds_realpower, bounds_imagpower, max_current
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



costs_generators, loads, shunts, bounds_voltage, bounds_realpower, bounds_imagpower, max_current =
get_features("D:\\repo\\data\\data_Matpower\\matpower\\case1888rte.m")
println(bounds_realpower["BUS_802"])
println(bounds_imagpower["BUS_802"])
# println(loads)
