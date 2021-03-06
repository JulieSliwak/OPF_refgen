export export_SDP

function export_SDP(sdp::SDPPrimal, path; indentedprint=true, renamemoments=true)

    # Build moment shortname dict if required
    momentdict = build_momentdict(sdp, renamemoments)

    # Collect all constraints keys
    ctr_keys = build_ctrkeysset(sdp)

    # Export blocks of constraints
    blocks_file = joinpath(path, "matrix.sdp")
    !isfile(blocks_file) || rm(blocks_file)

    fblocks = open(blocks_file, "a")
    print_blocksfile(fblocks, sdp.blocks, momentdict=momentdict, indentedprint=indentedprint)
    close(fblocks)

    # Export linear part of constraints
    lin_file = joinpath(path, "lin.sdp")
    !isfile(lin_file) || rm(lin_file)

    flin = open(lin_file, "a")
    print_linfile(flin, sdp.lin, sdp.linsym, momentdict=momentdict, indentedprint=indentedprint)
    close(flin)

    # Export constants of constraints
    cst_file = joinpath(path, "const.sdp")
    !isfile(cst_file) || rm(cst_file)

    fcst = open(cst_file, "a")
    print_cstfile(fcst, sdp.cst, momentdict=momentdict, ctr_keys=ctr_keys, indentedprint=indentedprint)
    close(fcst)


    # Export bloc types
    types_file = joinpath(path, "types.sdp")
    !isfile(types_file) || rm(types_file)

    ftypes = open(types_file, "a")
    print_typesfile(ftypes, sdp.block_to_vartype)
    close(ftypes)

    # Export moment dictionary name if required
    momname_file = joinpath(path, "names.sdp")
    !isfile(momname_file) || rm(momname_file)

    if renamemoments
        fname = open(momname_file, "a")
        print_namesfile(fname, momentdict)
        close(fname)
    end

    return nothing
end


"""
    momentdict = build_momentdict(sdp, renamemoments)

    Build a dict `momentdict::Dict{Exponent, String}` that stores short names for moments (Exponent) of `sdp`.
    NOTE: should be extended to Variable and (String, Variable) later...
"""
function build_momentdict(sdp, renamemoments::Bool)
    momentdict = Dict{Exponent, String}()

    n_moment=0
    n_matvar=0
    n_scalvar=0

    momentdict[Exponent()] = "1"
    ## Treating blocks
    for (moment, blockname, ??, ??) in keys(sdp.blocks)
        ??, ?? = moment.conj_part, moment.expl_part
        if renamemoments
            haskey(momentdict, ??) || (n_moment+=1; momentdict[??] = shortname_moment(n_moment))
            haskey(momentdict, ??) || (n_moment+=1; momentdict[??] = shortname_moment(n_moment))
            haskey(momentdict, ??) || (n_moment+=1; momentdict[??] = shortname_moment(n_moment))
            haskey(momentdict, ??) || (n_moment+=1; momentdict[??] = shortname_moment(n_moment))
        else
            haskey(momentdict, ??) || (momentdict[??] = format_string(??))
            haskey(momentdict, ??) || (momentdict[??] = format_string(??))
            haskey(momentdict, ??) || (momentdict[??] = format_string(??))
            haskey(momentdict, ??) || (momentdict[??] = format_string(??))
        end
    end

    for (moment, var) in keys(sdp.lin)
        ??, ?? = moment.conj_part, moment.expl_part
        if renamemoments
            haskey(momentdict, ??) || (n_moment+=1; momentdict[??] = shortname_moment(n_moment))
            haskey(momentdict, ??) || (n_moment+=1; momentdict[??] = shortname_moment(n_moment))
        else
            momentdict[??] = format_string(??)
            momentdict[??] = format_string(??)
        end
    end

    for (moment, blockname, var) in keys(sdp.linsym)
        ??, ?? = moment.conj_part, moment.expl_part
        if renamemoments
            haskey(momentdict, ??) || (n_moment+=1; momentdict[??] = shortname_moment(n_moment))
            haskey(momentdict, ??) || (n_moment+=1; momentdict[??] = shortname_moment(n_moment))
        else
            momentdict[??] = format_string(??)
            momentdict[??] = format_string(??)
        end
    end

    for moment in keys(sdp.cst)
        ??, ?? = moment.conj_part, moment.expl_part
        if renamemoments
            haskey(momentdict, ??) || (n_moment+=1; momentdict[??] = shortname_moment(n_moment))
            haskey(momentdict, ??) || (n_moment+=1; momentdict[??] = shortname_moment(n_moment))
        else
            momentdict[??] = format_string(??)
            momentdict[??] = format_string(??)
        end
    end

    return momentdict
end

function build_ctrkeysset(sdp::SDPPrimal{T}) where T
    ctr_keys = Set{Moment}()

    for ((moment, blockname, ??, ??), ??) in sdp.blocks
        push!(ctr_keys, moment)
    end
    for ((moment, var), ??) in sdp.lin
        push!(ctr_keys, moment)
    end
    for ((moment, blockname, var), ??) in sdp.linsym
        push!(ctr_keys, moment)
    end
    for (moment, ??) in sdp.cst
        push!(ctr_keys, moment)
    end
    return ctr_keys
end


function print_blocksfile(io::IO, sdpblocks::Dict{Tuple{Moment, String, Exponent, Exponent}, T};
                                                        momentdict::Dict{Exponent, String}=Dict{Exponent, String}(),
                                                        indentedprint=false,
                                                        print_header=true) where T
    if print_header
        println(io, "## Description of the matrices A_ji for the problem:")
        println(io, "##         max     ??? A_0i[k,l] ?? Zi[k,l] + ??? b_0[k] ?? x[k] + c_0")
        println(io, "##         s.t.    ??? A_ji[k,l] ?? Zi[k,l] + ??? b_j[k] ?? x[k] + c_j  ==  0")
        println(io, "## Constraints keys are j ??? (j_conj, j_expl, clique).")
        println(io, "## Objective keys are 0 ??? (1,1, *) for any *.")
        println(io, "#")
    end

    cstrlen?? = maximum(x->length(haskey(momentdict, x[1].conj_part) ? momentdict[x[1].conj_part] : string(x[1].conj_part)), keys(sdpblocks))
    cstrlen??= max(cstrlen??, length("#j_conj"))
    cstrlen?? = maximum(x->length(haskey(momentdict, x[1].expl_part) ? momentdict[x[1].expl_part] : string(x[1].expl_part)), keys(sdpblocks))
    cstrlen??= max(cstrlen??, length("j_expl"))
    cliquelen = maximum(x->length(x[1].clique), keys(sdpblocks))
    cliquelen= max(cliquelen, length("clique"))
    blocklen = maximum(x->length(x[2]), keys(sdpblocks))
    blocklen= max(blocklen, length("Zi"))
    rowlen = maximum(x->length(haskey(momentdict, x[3]) ? momentdict[x[3]] : string(x[3])), keys(sdpblocks))
    rowlen = max(rowlen, length("k"))
    collen = maximum(x->length(haskey(momentdict, x[4]) ? momentdict[x[4]] : string(x[4])), keys(sdpblocks))
    collen = max(collen, length("l"))

    print_string(io, "#j_conj", cstrlen??, indentedprint=indentedprint)
    print_string(io, "j_expl", cstrlen??, indentedprint=indentedprint)
    print_string(io, "clique", cliquelen, indentedprint=indentedprint)
    print_string(io, "Zi", blocklen, indentedprint=indentedprint)
    print_string(io, "k", rowlen, indentedprint=indentedprint)
    print_string(io, "l", collen, indentedprint=indentedprint)
    print_string(io, "Real(A_ij[k, l])", 23, indentedprint=indentedprint)
    print_string(io, "Imag(A_ij[k, l])", 23, indentedprint=indentedprint)
    println(io)

    for (moment, blockname, ??, ??) in sort(collect(keys(sdpblocks)))
        ?? = sdpblocks[(moment, blockname, ??, ??)]

        ??, ?? = moment.conj_part, moment.expl_part
        print_string(io, haskey(momentdict, ??) ? momentdict[??] : string(??), cstrlen??, indentedprint=indentedprint)
        print_string(io, haskey(momentdict, ??) ? momentdict[??] : string(??), cstrlen??, indentedprint=indentedprint)
        print_string(io, moment.clique, cliquelen, indentedprint=indentedprint)
        print_string(io, blockname, blocklen, indentedprint=indentedprint)
        print_string(io, haskey(momentdict, ??) ? momentdict[??] : string(??), rowlen, indentedprint=indentedprint)
        print_string(io, haskey(momentdict, ??) ? momentdict[??] : string(??), collen, indentedprint=indentedprint)
        @printf(io, "% .16e % .16e\n", real(??), imag(??))
    end
end


function print_linfile(io::IO, sdplin::Dict{Tuple{Moment, Exponent}, T}, sdplinsym::Dict{Tuple{Moment, String, Exponent}, T};
                                                                     momentdict::Dict{Exponent, String}=Dict{Exponent, String}(),
                                                                     indentedprint=false,
                                                                     print_header=true) where T
    if print_header
        println(io, "## Description of the vectors b_j for the problem:")
        println(io, "##         max     ??? A_0i[k,l] ?? Zi[k,l] + ??? b_0[k] ?? x[k] + c_0")
        println(io, "##         s.t.    ??? A_ji[k,l] ?? Zi[k,l] + ??? b_j[k] ?? x[k] + c_j  ==  0")
        println(io, "## Constraints keys are j ??? (j_conj, j_expl, clique).")
        println(io, "## Objective keys are 0 ??? (1,1, *) for any *.")
        println(io, "#")
    end

    cstrlen?? = 0
    cstrlen?? = 0
    cliquelen = 0
    varlen = varlensym = 0

    if length(union(keys(sdplin), keys(sdplinsym))) > 0
        cstrlen?? = maximum(x->length(haskey(momentdict, x[1].conj_part) ? momentdict[x[1].conj_part] : string(x[1].conj_part)), union(keys(sdplin), keys(sdplinsym)))
        cstrlen?? = maximum(x->length(haskey(momentdict, x[1].expl_part) ? momentdict[x[1].expl_part] : string(x[1].expl_part)), union(keys(sdplin), keys(sdplinsym)))
        cliquelen = maximum(x->length(x[1].clique), union(keys(sdplin), keys(sdplinsym)))
        varlen = length(sdplin)!=0 ? maximum(x->length(format_string(x[2])), keys(sdplin)) : 0
        varlensym = length(sdplinsym)!=0 ? maximum(x->length(format_string(x[3], x[2])), keys(sdplinsym)) : 0
    end

    cstrlen??= max(cstrlen??, length("#j_conj"))
    cstrlen??= max(cstrlen??, length("j_expl"))
    cliquelen= max(cliquelen, length("clique"))
    varlen = max(varlen, varlensym, length("x[k]"))

    print_string(io, "#j_conj", cstrlen??, indentedprint=indentedprint)
    print_string(io, "j_expl", cstrlen??, indentedprint=indentedprint)
    print_string(io, "clique", cliquelen, indentedprint=indentedprint)
    print_string(io, "x[k]", varlen, indentedprint=indentedprint)
    print_string(io, "Real(b_j[k])", 23, indentedprint=indentedprint)
    print_string(io, "Imag(b_j[k])", 23, indentedprint=indentedprint)
    println(io)

    if length(sdplin)!=0
        for (moment, var) in sort(collect(keys(sdplin)))
            ?? = sdplin[(moment, var)]

            ??, ?? = moment.conj_part, moment.expl_part
            print_string(io, haskey(momentdict, ??) ? momentdict[??] : string(??), cstrlen??, indentedprint=indentedprint)
            print_string(io, haskey(momentdict, ??) ? momentdict[??] : string(??), cstrlen??, indentedprint=indentedprint)
            print_string(io, moment.clique, cliquelen, indentedprint=indentedprint)
            print_string(io, format_string(var), varlen, indentedprint=indentedprint)
            @printf(io, "% .16e % .16e\n", real(??), imag(??))
        end
    end
    if length(sdplinsym) != 0
        for (moment, blockname, var) in sort(collect(keys(sdplinsym)))
            ?? = sdplinsym[(moment, blockname, var)]

            ??, ?? = moment.conj_part, moment.expl_part
            print_string(io, haskey(momentdict, ??) ? momentdict[??] : string(??), cstrlen??, indentedprint=indentedprint)
            print_string(io, haskey(momentdict, ??) ? momentdict[??] : string(??), cstrlen??, indentedprint=indentedprint)
            print_string(io, moment.clique, cliquelen, indentedprint=indentedprint)
            print_string(io, format_string(var, blockname), varlen, indentedprint=indentedprint)
            @printf(io, "% .16e % .16e\n", real(??), imag(??))
        end
    end
end


function print_cstfile(io::IO, sdpcst::Dict{Moment, T};
                                momentdict::Dict{Exponent, String}=Dict{Exponent, String}(),
                                ctr_keys::Set{Moment}=Set{Moment}(),
                                indentedprint=false,
                                print_header=true) where T
    if print_header
        println(io, "## Description of the scalars c_j for the problem:")
        println(io, "##         max     ??? A_0i[k,l] ?? Zi[k,l] + ??? b_0[k] ?? x[k] + c_0")
        println(io, "##         s.t.    ??? A_ji[k,l] ?? Zi[k,l] + ??? b_j[k] ?? x[k] + c_j  ==  0")
        println(io, "## Constraints keys are j ??? (j_conj, j_expl, clique).")
        println(io, "## Objective keys are 0 ??? (1,1, *) for any *.")
        println(io, "#")
    end

    cstrlen?? = maximum(x->length(haskey(momentdict, x.conj_part) ? momentdict[x.conj_part] : string(x.conj_part)), union(ctr_keys, keys(sdpcst)))
    cstrlen??= max(cstrlen??, length("#j_conj"))
    cstrlen?? = maximum(x->length(haskey(momentdict, x.expl_part) ? momentdict[x.expl_part] : string(x.expl_part)), union(ctr_keys, keys(sdpcst)))
    cstrlen??= max(cstrlen??, length("j_expl"))
    cliquelen = maximum(x->length(x.clique), union(ctr_keys, keys(sdpcst)))
    cliquelen= max(cliquelen, length("clique"))

    print_string(io, "#j_conj", cstrlen??, indentedprint=indentedprint)
    print_string(io, "j_expl", cstrlen??, indentedprint=indentedprint)
    print_string(io, "clique", cliquelen, indentedprint=indentedprint)
    print_string(io, "Real(c_j)", 23, indentedprint=indentedprint)
    print_string(io, "Imag(c_j)", 23, indentedprint=indentedprint)
    println(io)

    for moment in union(ctr_keys, keys(sdpcst))
    # for moment in sort(union(ctr_keys, keys(sdpcst)))
        ??, ?? = moment.conj_part, moment.expl_part
        ?? = haskey(sdpcst, moment) ? sdpcst[moment] : 0
        print_string(io, haskey(momentdict, ??) ? momentdict[??] : string(??), cstrlen??, indentedprint=indentedprint)
        print_string(io, haskey(momentdict, ??) ? momentdict[??] : string(??), cstrlen??, indentedprint=indentedprint)
        print_string(io, moment.clique, cliquelen, indentedprint=indentedprint)
        @printf(io, "% .16e % .16e\n", real(??), imag(??))
    end
end


function print_typesfile(io::IO, block_to_vartype::Dict{String, Symbol})
    println(io, "## Description of the matrix and scalar variables Zi and x[k] for the problem:")
    println(io, "##         max     ??? A_0i[k,l] ?? Zi[k,l] + ??? b_0[k] ?? x[k] + c_0")
    println(io, "##         s.t.    ??? A_ji[k,l] ?? Zi[k,l] + ??? b_j[k] ?? x[k] + c_j  ==  0")
    println(io, "## Matrix variable types Zi are \"SDPC\" or \"SDP\".")
    println(io, "## By default, scalar variables x[k] are assumed to be real, free.")
    println(io, "#")

    cstrlen = maximum(x->length(x), keys(block_to_vartype))
    cstrlen = max(cstrlen, length("#Zi"))

    print_string(io, "#Zi", cstrlen, alignright=false); println(io, " type")

    for blockname in sort(collect(keys(block_to_vartype)))
        vartype = block_to_vartype[blockname]
        if vartype in Set([:SDP, :SDPC])
            print_string(io, blockname, cstrlen)
            println(io, " $(string(vartype))")
        end
        # warn(LOGGER, "Ignoring variable $blockname of type $vartype") ## NOTE: better logging system.
    end
end


function print_namesfile(io::IO, momentdict::Dict{Exponent, String})
    println(io, "## Description of the scalars c_j for the problem:")
    println(io, "##         max     ??? A_0i[k,l] ?? Zi[k,l] + ??? b_0[k] ?? x[k] + c_0")
    println(io, "##         s.t.    ??? A_ji[k,l] ?? Zi[k,l] + ??? b_j[k] ?? x[k] + c_j  ==  0")
    println(io, "## Constraints keys are j ??? (j_conj, j_expl, clique).")
    println(io, "## Objective keys are 0 ??? (1,1, *) for any *.")
    println(io, "#")
    println(io, "#shortname  Explicit_name")

    for ?? in sort(collect(keys(momentdict)))
        shortname = momentdict[??]
        println(io, "$shortname $(format_string(??))")
    end
end
