using LightGraphs, MetaGraphs

function read_network_graph(instance_path::String)

    data = load_matpower(instance_path)

    ## Bus load and shunt information
    i_debut, i_fin = find_numarray(1, data)
    checkfor(data, i_debut-1, "mpc.bus")
    nb_bus = i_fin-i_debut+1
    index_bus = Dict( data[i+i_debut-1,1] => i for i=1:nb_bus)
    if index_bus != Dict( i => i for i=1:nb_bus)
        println("!!! Specific numerotation of buses !!! \n")
    end
    ## Bus generator information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.gen")

    #initialize network graph G
    G = MetaGraph(nb_bus)
    for i in 1:nb_bus
        set_props!(G, i, Dict(:name => "node$i"))
    end

    ## Link information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.branch")
    for i=i_debut:i_fin
      if data[i, 11] == 0
        #@warn("link $(data[i,1])⟶$(data[i,2]) breaker out of service !")
      else
          orig = index_bus[data[i,1]]
          dest = index_bus[data[i,2]]
          LightGraphs.add_edge!(G, orig, dest)
      end
    end
    return G
end


function read_network_graph_real_numbers(instance_path::String)

    data = load_matpower(instance_path)

    ## Bus load and shunt information
    i_debut, i_fin = find_numarray(1, data)
    checkfor(data, i_debut-1, "mpc.bus")
    nb_bus = i_fin-i_debut+1
    index_bus = Dict( data[i+i_debut-1,1] => i for i=1:nb_bus)
    if index_bus != Dict( i => i for i=1:nb_bus)
        println("!!! Specific numerotation of buses !!! \n")
    end
    ## Bus generator information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.gen")

    #initialize network graph G
    G = MetaGraph(2*nb_bus)

    #real parts
    for i in 1:nb_bus
        set_props!(G, i, Dict(:name => "node$(i)"))  #bus i => node i for iRe, node nb_bus+i for iIm
    end

    #imaginary parts
    for i in nb_bus+1:2*nb_bus
        set_props!(G, i, Dict(:name => "node$(i)"))
    end

    ## Link information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.branch")
    for i=i_debut:i_fin
      if data[i, 11] == 0
        #@warn("link $(data[i,1])⟶$(data[i,2]) breaker out of service !")
      else
          orig = index_bus[data[i,1]]
          dest = index_bus[data[i,2]]
          node_orig_Re = orig
          node_orig_Im = nb_bus + orig
          node_dest_Re = dest
          node_dest_Im = nb_bus + dest
          add_edge!(G, node_orig_Re, node_dest_Re) #1 edge in complex = 4 edges in real
          add_edge!(G, node_orig_Re, node_dest_Im) #1 edge in complex = 4 edges in real
          add_edge!(G, node_orig_Im, node_dest_Re) #1 edge in complex = 4 edges in real
          add_edge!(G, node_orig_Im, node_dest_Im) #1 edge in complex = 4 edges in real
      end
    end
    return G
end


using LinearAlgebra, SparseArrays, SuiteSparse

function read_sparsity_pattern(instance_path::String)

    data = load_matpower(instance_path)

    ## Bus load and shunt information
    i_debut, i_fin = find_numarray(1, data)
    checkfor(data, i_debut-1, "mpc.bus")
    nb_bus = i_fin-i_debut+1
    index_bus = Dict( data[i+i_debut-1,1] => i for i=1:nb_bus)

    ## Bus generator information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.gen")

    #initialize network graph G
    sp = spzeros(nb_bus,nb_bus)

    ## Link information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.branch")
    for i=i_debut:i_fin
      if data[i, 11] == 0
        #@warn("link $(data[i,1])⟶$(data[i,2]) breaker out of service !")
      else
          orig = index_bus[data[i,1]]
          dest = index_bus[data[i,2]]
          sp[orig,dest] = 1
      end
    end

    sp_sym = sp + sp'

    diag = zeros(nb_bus)

    for i=1:nb_bus
        sum_col = sum(sp_sym[i,:])
        diag[i] = Int(sum_col + 1)
    end
    return sp_sym + Diagonal(diag)
    # return sp+sp'+nb_bus*sparse(I, nb_bus, nb_bus)
end


function read_sparsity_pattern_real(instance_path::String)
    data = load_matpower(instance_path)
    ## Bus load and shunt information
    i_debut, i_fin = find_numarray(1, data)
    checkfor(data, i_debut-1, "mpc.bus")
    nb_bus = i_fin-i_debut+1
    index_bus = Dict( data[i+i_debut-1,1] => i for i=1:nb_bus)

    ## Bus generator information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.gen")

    sp = spzeros(2*nb_bus,2*nb_bus)

    ## Link information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.branch")
    for i=i_debut:i_fin
      if data[i, 11] == 0
        #@warn("link $(data[i,1])⟶$(data[i,2]) breaker out of service !")
      else
          orig = index_bus[data[i,1]]
          dest = index_bus[data[i,2]]
          node_orig_Re = orig
          node_orig_Im = nb_bus + orig
          node_dest_Re = dest
          node_dest_Im = nb_bus + dest               #bus i => node i for iRe, node nb_bus+i for iIm
          sp[node_orig_Re, node_dest_Re] = 1
          sp[node_orig_Re, node_dest_Im] = 1
          sp[node_orig_Im, node_dest_Re] = 1
          sp[node_orig_Im, node_dest_Im] = 1
      end
    end
    sp_sym = sp + sp'
    diag = zeros(2*nb_bus)
    for i=1:2*nb_bus
        sum_col = sum(sp_sym[i,:])
        diag[i] = Int(sum_col + 1)
    end
    return sp_sym + Diagonal(diag)
end


function read_correlative_sparsity_graph(instance_path::String)

    data = load_matpower(instance_path)

    ## Bus load and shunt information
    i_debut, i_fin = find_numarray(1, data)
    checkfor(data, i_debut-1, "mpc.bus")
    nb_bus = i_fin-i_debut+1
    index_bus = Dict( data[i+i_debut-1,1] => i for i=1:nb_bus)
    if index_bus != Dict( i => i for i=1:nb_bus)
        println("!!! Specific numerotation of buses !!! \n")
    end
    ## Bus generator information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.gen")

    #initialize network graph G
    G = MetaGraph(nb_bus)
    for i in 1:nb_bus
        set_props!(G, i, Dict(:name => "node$i"))
    end

    ## Link information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.branch")
    for i=i_debut:i_fin
      if data[i, 11] == 0
        #@warn("link $(data[i,1])⟶$(data[i,2]) breaker out of service !")
      else
          orig = index_bus[data[i,1]]
          dest = index_bus[data[i,2]]
          LightGraphs.add_edge!(G, orig, dest)
      end
    end

    C =  MetaGraph(nb_bus)
    for i in 1:nb_bus
        set_props!(C, i, Dict(:name => "node$i"))
    end

    for i in 1:nb_bus
        neighbors_i = neighbors(G,i)
        for u in union(i, neighbors_i)
            for v in union(i, neighbors_i)
                if !has_edge(C, u, v) && u != v
                    LightGraphs.add_edge!(C, u, v)
                end
            end
        end
    end
    return C
end

function read_correlative_sparsity_graph_real(instance_path::String)

    data = load_matpower(instance_path)

    ## Bus load and shunt information
    i_debut, i_fin = find_numarray(1, data)
    checkfor(data, i_debut-1, "mpc.bus")
    nb_bus = i_fin-i_debut+1
    index_bus = Dict( data[i+i_debut-1,1] => i for i=1:nb_bus)
    if index_bus != Dict( i => i for i=1:nb_bus)
        println("!!! Specific numerotation of buses !!! \n")
    end
    ## Bus generator information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.gen")

    #initialize network graph G
    G = MetaGraph(2*nb_bus)
    for i in 1:nb_bus
        set_props!(G, i, Dict(:name => "node$i"))
    end

    ## Link information
    i_debut, i_fin = find_numarray(i_fin+1, data)
    checkfor(data, i_debut-1, "mpc.branch")
    for i=i_debut:i_fin
      if data[i, 11] == 0
        #@warn("link $(data[i,1])⟶$(data[i,2]) breaker out of service !")
      else
          orig = index_bus[data[i,1]]
          dest = index_bus[data[i,2]]
          node_orig_Re = orig
          node_orig_Im = nb_bus + orig
          node_dest_Re = dest
          node_dest_Im = nb_bus + dest               #bus i => node i for iRe, node nb_bus+i for iIm
          LightGraphs.add_edge!(G, node_orig_Re, node_dest_Re)
          LightGraphs.add_edge!(G, node_orig_Re, node_dest_Im)
          LightGraphs.add_edge!(G, node_orig_Im, node_dest_Re)
          LightGraphs.add_edge!(G, node_orig_Im, node_dest_Im)
      end
    end

    C =  MetaGraph(2*nb_bus)
    for i in 1:2*nb_bus
        set_props!(C, i, Dict(:name => "node$i"))
    end

    for i in 1:2*nb_bus
        neighbors_i = neighbors(G,i)
        for u in union(i, neighbors_i)
            for v in union(i, neighbors_i)
                if !has_edge(C, u, v) && u != v
                    LightGraphs.add_edge!(C, u, v)
                end
            end
        end
    end
    return C
end
