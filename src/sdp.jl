function certify_gap(N::Int, H::ncpoly, gamma, d::Int; QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    ptsupp = get_wbasis(N, min(N, 2d))
    ptsupp = ptsupp[2:end]
    sort!(ptsupp, lt=isless_td)
    tbasis,wbasis,sbasis = get_tbasis(N, d, ptsupp)
    lb = length(tbasis)
    gbasis = get_wbasis(N, N-1)
    lgb = length(gbasis)
    if QUIET == false
        println("The block sizes are [$lb, $lgb].")
    end
    tsupp = Vector{Int}[]
    for i = 1:length(tbasis), j = i:length(tbasis)
        @inbounds bi1 = reduce!([wbasis[tbasis[i][1]]; wbasis[tbasis[j][1]]])[1]
        @inbounds bi2 = sort([sbasis[tbasis[i][2]]; sbasis[tbasis[j][2]]])
        bi = state_reduce(bi1, bi2, ptsupp)
        push!(tsupp, bi)
    end
    sort!(tsupp)
    unique!(tsupp)
    ltsupp = length(tsupp)
    if QUIET == false
        println("There are $ltsupp affine constraints.")
        println("Assembling the SDP...")
    end
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:ltsupp]
    @inbounds pos = @variable(model, [1:2*lb, 1:2*lb], PSD)
    for i = 1:lb, j = i:lb
        @inbounds bi1,c = reduce!([wbasis[tbasis[i][1]]; wbasis[tbasis[j][1]]])
        @inbounds bi2 = sort([sbasis[tbasis[i][2]]; sbasis[tbasis[j][2]]])
        bi = state_reduce(bi1, bi2, ptsupp)
        Locb = bfind(tsupp, ltsupp, bi)
        pp1 = pos[i, j] + pos[i+lb, j+lb]
        pp2 = pos[i+lb, j] - pos[j+lb, i]
        if i == j
            @inbounds add_to_expression!(cons[Locb], real(c), pp1)
        else
            @inbounds add_to_expression!(cons[Locb], 2*real(c), pp1)
            if imag(c) != 0
                @inbounds add_to_expression!(cons[Locb], -2*imag(c), pp2)
            end
        end
    end
    @inbounds gpos = @variable(model, [1:2*lgb, 1:2*lgb], PSD)
    for i = 1:lgb, j = i:lgb
        pp1 = gpos[i, j] + gpos[i+lgb, j+lgb]
        pp2 = gpos[i+lgb, j] - gpos[j+lgb, i]
        for k = 1:length(H.supp)
            @inbounds bi,c = reduce!([gbasis[i]; H.supp[k]; gbasis[j]])
            if isempty(bi)
                Locb = 1
            else
                loc = bfind(ptsupp, length(ptsupp), bi, lt=isless_td)
                Locb = bfind(tsupp, ltsupp, [loc])
            end
            if i == j
                @inbounds add_to_expression!(cons[Locb], real(c), pp1)
            else
                @inbounds add_to_expression!(cons[Locb], 2*real(c), pp1)
                if imag(c) != 0
                    @inbounds add_to_expression!(cons[Locb], -2*imag(c), pp2)
                end
            end
            @inbounds bi,c = reduce!([gbasis[i]; gbasis[j]; H.supp[k]])
            if isempty(bi)
                Locb = 1
            else
                loc = bfind(ptsupp, length(ptsupp), bi, lt=isless_td)
                Locb = bfind(tsupp, ltsupp, [loc])
            end
            if i == j
                @inbounds add_to_expression!(cons[Locb], -real(c), pp1)
            else
                @inbounds add_to_expression!(cons[Locb], -2*real(c), pp1)
                if imag(c) != 0
                    @inbounds add_to_expression!(cons[Locb], 2*imag(c), pp2)
                end
            end
        end
        @inbounds bi,c = reduce!([gbasis[i]; gbasis[j]])
        if isempty(bi)
            Locb = 1
        else
            loc = bfind(ptsupp, length(ptsupp), bi, lt=isless_td)
            Locb = bfind(tsupp, ltsupp, [loc])
        end
        if i == j
            @inbounds add_to_expression!(cons[Locb], -gamma*real(c), pp1)
        else
            @inbounds add_to_expression!(cons[Locb], -2*gamma*real(c), pp1)
            if imag(c) != 0
                @inbounds add_to_expression!(cons[Locb], 2*gamma*imag(c), pp2)
            end
        end
        if isempty(gbasis[i])
            loc1 = Int[]
        else
            loc1 = bfind(ptsupp, length(ptsupp), gbasis[i], lt=isless_td)
        end
        if isempty(gbasis[j])
            loc2 = Int[]
        else
            loc2 = bfind(ptsupp, length(ptsupp), gbasis[j], lt=isless_td)
        end
        Locb = bfind(tsupp, ltsupp, sort([loc1; loc2]))
        if i == j
            @inbounds add_to_expression!(cons[Locb], gamma, pp1)
        else
            @inbounds add_to_expression!(cons[Locb], 2*gamma, pp1)
        end
    end
    @variable(model, lower)
    cons[1] += lower
    @constraint(model, cons .== 0)
    @objective(model, Max, lower)
    if QUIET == false
        println("Solving the SDP...")
    end
    time = @elapsed begin
        optimize!(model)
    end
    if QUIET == false
        println("SDP solving time: $time seconds.")
    end
    status = termination_status(model)
    @show status
    return status
end
