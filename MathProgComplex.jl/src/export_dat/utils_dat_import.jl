export import_from_dat

"""
    pb, point = import_from_dat(instancepath::String, precondcstrspath::String)

Build the polynomial optimization problem described by the `instancepath` file,
with the dispensory preconditionning descritpion `precondcstrspath` along with
the initial point possibly provided in the file (defaults value is null).
"""
function import_from_dat(instancepath::String, precondfilename::String="")
    point = Point()
    variables = SortedDict{String, Variable}()
    exponents = SortedDict{String, Exponent}()
    pb = Problem()

    instance_str = open(instancepath)
    l = jump_comments!(instance_str)


    ## Collect and define variables
    # line = matchall(r"\S+", l)
    line = collect(( m.match for m=eachmatch(r"\S+", l) ))
    while line[1] == "VAR_TYPE" && !eof(instance_str)
        if line[2] == "REAL"
            var = Variable(line[3], Real)
        elseif line[2] == "BOOL"
            var = Variable(line[3], Bool)
        elseif line[2] == "C"
            var = Variable(line[3], Complex)
        else
            error(LOGGER, "import_to_dat(): Unknown variable type $(line[2]) for variable $(line[3]).")
        end

        variables[line[3]] = var
        add_variable!(pb, var)

        setindex!(point, convert(var.kind, parse.(Float64, line[5]) + im*parse.(Float64,line[6])), var)

        ## Mark where objective definition begins
        mark(instance_str)

        l = readline(instance_str)
        line = collect(( m.match for m=eachmatch(r"\S+", l) ))
    end


    ## Move forward to the monomial definition section
    while line[1] != "MONO_DEF" && !eof(instance_str)
        l = readline(instance_str)
        line = collect(( m.match for m=eachmatch(r"\S+", l) ))
    end

    ## Build all monomials
    while line[1] == "MONO_DEF" && !eof(instance_str)
        exponame = line[2]
        var = variables[line[3]]
        if !haskey(exponents, exponame)
            exponents[exponame] = Exponent()
        end
        exponents[exponame] = product(exponents[exponame], Exponent(SortedDict(var=>Degree(parse.(Float64,line[5]), parse.(Float64,line[6])))))
        l = readline(instance_str)
        line = collect(( m.match for m=eachmatch(r"\S+", l) ))
    end

    ## Add final monomial
    if line[1] == "MONO_DEF"
        exponame = line[2]
        var = variables[line[3]]
        if !haskey(exponents, exponame)
            exponents[exponame] = Exponent()
        end
        exponents[exponame] = product(exponents[exponame], Exponent(SortedDict(var=>Degree(parse.(Float64,line[5]), parse.(Float64,line[6])))))
    end

    ## Reset stream to the objective defintion
    reset(instance_str)
    l = readline(instance_str)
    line = collect(( m.match for m=eachmatch(r"\S+", l) ))

    ## Build polynomial objective
    p = Polynomial()
    while line[2] == "OBJ"
        ?? = parse_??(line[5], line[6])
        var1, var2 = line[3:4]
        if line[1] == "MONO"
            add!(p, ?? * exponents[var1])
        else
            quad_expo = Exponent()
            (var1!="NONE") && product!(quad_expo, conj(variables[var1]))
            (var2!="NONE") && product!(quad_expo, variables[var2])
            add!(p, Polynomial(SortedDict{Exponent, Number}(quad_expo=>??)))
        end
        l = readline(instance_str)
        line = collect(( m.match for m=eachmatch(r"\S+", l) ))
    end
    set_objective!(pb, p)

    ## Build constraints
    next_state = :AssembleCtr
    cur_ctr = line[2]
    var1, var2 = line[3:4]

    lb = -Inf
    ub = +Inf
    p = Polynomial()
    ?? = parse_??(line[5], line[6])
    while next_state!=:ReadLine || (!eof(instance_str) && line[1] != "MONO_DEF")

        state = next_state
        if state == :ReadLine
            ## Readline
            l = readline(instance_str)
            line = collect(( m.match for m=eachmatch(r"\S+", l) ))

            var1, var2 = line[3:4]
            ?? = parse_??(line[5], line[6])

            ## Determine whether current ctr should be completed or new ctr should be set
            if line[2] == cur_ctr
                next_state = :AssembleCtr
            else
                next_state = :SaveCtr
            end

            ## If monomial difnition section reached, exit loop
            if line[1] == "MONO_DEF"
                next_state = :ReadLine
            end

        elseif state == :AssembleCtr
            if line[1] == "MONO"
                add!(p, ?? * exponents[var1])
            elseif line[1] ??? SortedSet(["CONST", "LIN", "QUAD"])
                quad_expo = Exponent()
                (var1!="NONE") && product!(quad_expo, conj(variables[var1]))
                (var2!="NONE") && product!(quad_expo, variables[var2])
                add!(p, Polynomial(SortedDict{Exponent, Number}(quad_expo=>??)))
            elseif line[1] == "UB"
                ub = ??
            elseif line[1] == "LB"
                lb = ??
            else
                error(LOGGER, "import_from_dat(): Unknown variable type $(line[1]) for constraint $(line[2]).")
            end

            next_state = :ReadLine

        elseif state == :SaveCtr
            add_constraint!(pb, String(cur_ctr), lb << p << ub)

            lb = -Inf
            ub = +Inf
            p = Polynomial()
            cur_ctr = line[2]
            next_state = :AssembleCtr
        else
            error(LOGGER, "import_from_dat(): Unknown state $state")
        end
    end

    # Add final constraint
    add_constraint!(pb, String(cur_ctr), lb << p << ub)

    close(instance_str)

    ## Set preconditioning flag
    if precondfilename != ""
        precond_str = open(joinpath(splitdir(instancepath)[1], precondfilename))

        l = jump_comments!(precond_str)
        line = collect(( m.match for m=eachmatch(r"\S+", l) ))
        while !eof(precond_str)
            if line[2] == "SQRT"
                pb.constraints[line[1]].precond = :sqrt
            else
                error(LOGGER, "import_from_dat(): Unknown preconditioning $(line[2]) for constraint $(line[1]).")
            end
            l = readline(precond_str)
            line = collect(( m.match for m=eachmatch(r"\S+", l) ))
        end
        close(precond_str)
    end

    return pb, point
end


"""
    l = jump_comments!(io::IOStream)

Jump comments, i.e. lines starting by '#', in the `io` stream, and return the
first non commented line.
"""
function jump_comments!(io::IOStream)
    l = ""
    while !eof(io)
        l = readline(io)
        occursin(r"\s*#", l) || break
    end
    return l
end

function parse_??(realpart::T, imagpart::T) where T <: Union{String, SubString}
    ?? = 0
    if realpart == "Inf"
        ?? += Inf
    elseif realpart == "-Inf"
        ?? += -Inf
    else
        ?? += parse.(Float64,realpart)
    end

    if imagpart == "Inf"
        ?? += Inf*im
    elseif imagpart == "-Inf"
        ?? += -Inf*im
    else
        ?? += parse.(Float64,imagpart)*im
    end
    return ??
end
