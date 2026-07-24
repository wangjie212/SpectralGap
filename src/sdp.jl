function certify_Ising_gap(N::Int, H::ncpoly, gamma, d::Int; lso=6, QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    basis = [get_basis(N, d, label=i) for i in [1,2]]
    lb = length.(basis)
    gbasis = [get_bulkbasis(N, d-1, label=i) for i in [1,2]] 
    lgb = length.(gbasis)
    bs = [lb; lgb]
    if QUIET == false
        println("The block sizes are $bs.")
    end
    tsupp = Vector{Vector{Int}}[]
    for i = 1:length(basis), j = 1:lb[i], k = j:lb[i]
        @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]], N)
        if c != 0
            if isempty(bi)
                push!(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
            else
                push!(tsupp, sort([basis[i][j][2]; basis[i][k][2]; [bi]]))
            end
        end
    end
    for l = 1:length(gbasis), i = 1:lgb[l], j = i:lgb[l]
        bis = PSDstate_entry(gbasis[l][i][1], gbasis[l][j][1], H, N)[1]
        for bi in bis
            if isempty(bi)
                push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
            else
                push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
            end
        end
        if !isz(gbasis[l][i][1], model="Ising") && !isz(gbasis[l][j][1], model="Ising")
            temp = [reduce_mirror(gbasis[l][i][1], N), reduce_mirror(gbasis[l][j][1], N)]
            push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; temp]))
        end
    end
    # for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3
    #     ind = [i,j,k,l,s,t,u]
    #     inx = ind .!= 0
    #     bi = 3*(Vector(1:7)[inx] .- 1) + ind[inx]
    #     if !isz(bi)
    #         push!(tsupp, [reduce_mirror(bi, N)])
    #     end
    # end
    sort!(tsupp)
    unique!(tsupp)
    if QUIET == false
        println("There are $(length(tsupp)) affine constraints.")
        println("Assembling the SDP...")
    end
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:length(tsupp)]
    pos = Vector{Symmetric{VariableRef}}(undef, length(basis))
    for i = 1:length(basis)
        pos[i] = @variable(model, [1:lb[i], 1:lb[i]], PSD)
        for j = 1:lb[i], k = j:lb[i]
            @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]], N, realify=true)
            if c != 0
                if isempty(bi)
                    Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
                else
                    Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]; [bi]]))
                end
                if j == k
                    @inbounds add_to_expression!(cons[Locb], c, pos[i][j,k])
                else
                    @inbounds add_to_expression!(cons[Locb], 2c, pos[i][j,k])
                end
            end
        end
    end
    gpos = Vector{Symmetric{VariableRef}}(undef, length(gbasis))
    for l = 1:length(gbasis)
        gpos[l] = @variable(model, [1:lgb[l], 1:lgb[l]], PSD)
        for i = 1:lgb[l], j = i:lgb[l]
            bis,coes = PSDstate_entry(gbasis[l][i][1], gbasis[l][j][1], H, N, realify=true)
            for (k, bi) in enumerate(bis)
                if isempty(bi)
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
                else
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
                end
                if i == j
                    @inbounds add_to_expression!(cons[Locb], coes[k], gpos[l][i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], 2*coes[k], gpos[l][i, j])
                end
            end
            @inbounds bi,c = reduce!([gbasis[l][i][1]; gbasis[l][j][1]], N, realify=true)
            if c != 0
                if isempty(bi)
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
                else
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
                end
                if i == j
                    @inbounds add_to_expression!(cons[Locb], -c*gamma, gpos[l][i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], -2c*gamma, gpos[l][i, j])
                end
            end
            if !isz(gbasis[l][i][1], model="Ising") && !isz(gbasis[l][j][1], model="Ising")
                temp = [reduce_mirror(gbasis[l][i][1], N), reduce_mirror(gbasis[l][j][1], N)]
                Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; temp]))
                if i == j
                    @inbounds add_to_expression!(cons[Locb], gamma, gpos[l][i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], 2*gamma, gpos[l][i, j])
                end
            end
        end
    end
    mons = generate_mons(N, lso)
    mons = filter_mons(mons, tsupp, H, N)
    for mon in mons
        fr = @variable(model)
        for (i, item) in enumerate(H.supp)
            bi,c = reduce!([item; mon], N)
            if imag(c) != 0
                Locb = bfind(tsupp, [bi])
                add_to_expression!(cons[Locb], H.coe[i]*imag(c), fr)
            end
        end
    end
    # posepsd7!(model, cons, tsupp, N)
    @variable(model, λ)
    cons[1] += λ
    # for i = 1:ceil(Int, N/2)
    #     μ = @variable(model, lower_bound = 0)
    #     add_to_expression!(cons[1], μ)
    #     Locb = bfind(tsupp, [[3i-2], [3i-2]])
    #     add_to_expression!(cons[Locb], -μ)
    # end
    @objective(model, Max, λ)
    @constraint(model, cons .== 0)
    # @constraint(model, con, cons==zeros(length(cons)))
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
    # dual_var = -dual(con)
    # mmat = Vector{Matrix{ComplexF16}}(undef, length(basis))
    # for i = 1:length(basis)
    #     mmat[i] = zeros(ComplexF16, lb[i], lb[i])
    #     for j = 1:lb[i], k = 1:lb[i]
    #         @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]], N)
    #         if c != 0
    #             if isempty(bi)
    #                 Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
    #             else
    #                 Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]; [bi]]))
    #             end
    #             mmat[i][j,k] = c*dual_var[Locb]
    #         end
    #     end
    # end
    flag = status == MathOptInterface.OPTIMAL ? 1 : 0
    return flag
end

function certify_Ising_gap_nosignsymmetry(N::Int, H::ncpoly, gamma, d::Int; lso=6, QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    basis = [get_basis(N, d, label=1)]
    lb = length.(basis)
    gbasis = [[get_bulkbasis(N, d-1, label=1); get_bulkbasis(N, d-1, label=2)]]
    lgb = length.(gbasis)
    bs = [lb; lgb]
    if QUIET == false
        println("The block sizes are $bs.")
    end
    tsupp = Vector{Vector{Int}}[]
    for i = 1:length(basis), j = 1:lb[i], k = j:lb[i]
        @inbounds bi = reduce!([basis[i][j][1]; basis[i][k][1]], N, identify_zeros=false)[1]
        if isempty(bi)
            push!(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
        else
            push!(tsupp, sort([basis[i][j][2]; basis[i][k][2]; [bi]]))
        end
        temp = [reduce_mirror(basis[i][j][1], N), reduce_mirror(basis[i][k][1], N)]
        push!(tsupp, sort([basis[i][j][2]; basis[i][k][2]; temp]))
    end
    for l = 1:length(gbasis), i = 1:lgb[l], j = i:lgb[l]
        bis = PSDstate_entry(gbasis[l][i][1], gbasis[l][j][1], H, N, identify_zeros=false)[1]
        for bi in bis
            if isempty(bi)
                push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
            else
                push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
            end
        end
        temp = [reduce_mirror(gbasis[l][i][1], N), reduce_mirror(gbasis[l][j][1], N)]
        push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; temp]))
    end
    # for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3
    #     ind = [i,j,k,l,s,t,u]
    #     inx = ind .!= 0
    #     bi = 3*(Vector(1:7)[inx] .- 1) + ind[inx]
    #     push!(tsupp, [reduce_mirror(bi, N)])
    # end
    for i = 1:N, j = 1:3
        push!(tsupp, [[3(i-1)+j], [3(i-1)+j]])
    end
    sort!(tsupp)
    unique!(tsupp)
    if QUIET == false
        println("There are $(length(tsupp)) affine constraints.")
        println("Assembling the SDP...")
    end
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:length(tsupp)]
    pos = Vector{Symmetric{VariableRef}}(undef, length(basis))
    for i = 1:length(basis)
        pos[i] = @variable(model, [1:2*lb[i], 1:2*lb[i]], PSD)
        for j = 1:lb[i], k = j:lb[i]
            @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]], N, identify_zeros=false)
            if isempty(bi)
                Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
            else
                Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]; [bi]]))
            end
            if j == k
                @inbounds add_to_expression!(cons[Locb], c, pos[i][j,k]+pos[i][j+lb[i],k+lb[i]])
            elseif isreal(c)
                @inbounds add_to_expression!(cons[Locb], 2c, pos[i][j,k]+pos[i][j+lb[i],k+lb[i]])
            else
                @inbounds add_to_expression!(cons[Locb], -2*imag(c), pos[i][j,k+lb[i]]-pos[i][k,j+lb[i]])
            end
        end
    end
    npos = Vector{Symmetric{VariableRef}}(undef, length(basis))
    for l = 1:length(basis)
        npos[l] = @variable(model, [1:2*lb[l], 1:2*lb[l]], PSD)
        for i = 1:lb[l], j = i:lb[l]
            @inbounds bi,c = reduce!([basis[l][i][1]; basis[l][j][1]], N, identify_zeros=false)
            if isempty(bi)
                Locb = bfind(tsupp, sort([basis[l][i][2]; basis[l][j][2]]))
            else
                Locb = bfind(tsupp, sort([basis[l][i][2]; basis[l][j][2]; [bi]]))
            end
            if i == j
                @inbounds add_to_expression!(cons[Locb], c, npos[l][i, j]+npos[l][i+lb[l], j+lb[l]])
            elseif isreal(c)
                @inbounds add_to_expression!(cons[Locb], 2c, npos[l][i, j]+npos[l][i+lb[l], j+lb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], -2*imag(c), npos[l][i,j+lb[l]]-npos[l][j,i+lb[l]])
            end
            temp = [reduce_mirror(basis[l][i][1], N), reduce_mirror(basis[l][j][1], N)]
            Locb = bfind(tsupp, sort([basis[l][i][2]; basis[l][j][2]; temp]))
            if i == j
                @inbounds add_to_expression!(cons[Locb], -1, npos[l][i, j]+npos[l][i+lb[l], j+lb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], -2, npos[l][i, j]+npos[l][i+lb[l], j+lb[l]])
            end
        end
    end
    gpos = Vector{Symmetric{VariableRef}}(undef, length(gbasis))
    for l = 1:length(gbasis)
        gpos[l] = @variable(model, [1:2*lgb[l], 1:2*lgb[l]], PSD)
        for i = 1:lgb[l], j = i:lgb[l]
            bis,coes = PSDstate_entry(gbasis[l][i][1], gbasis[l][j][1], H, N, identify_zeros=false)
            for (k, bi) in enumerate(bis)
                if isempty(bi)
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
                else
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
                end
                if i == j
                    @inbounds add_to_expression!(cons[Locb], coes[k], gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
                elseif isreal(coes[k])
                    @inbounds add_to_expression!(cons[Locb], 2*coes[k], gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
                else
                    @inbounds add_to_expression!(cons[Locb], -2*imag(coes[k]), gpos[l][i, j+lgb[l]]-gpos[l][j, i+lgb[l]])
                end
            end
            @inbounds bi,c = reduce!([gbasis[l][i][1]; gbasis[l][j][1]], N, identify_zeros=false)
            if isempty(bi)
                Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
            else
                Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
            end
            if i == j
                @inbounds add_to_expression!(cons[Locb], -c*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            elseif isreal(c)
                @inbounds add_to_expression!(cons[Locb], -2c*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], 2*imag(c)*gamma, gpos[l][i, j+lgb[l]]-gpos[l][j, i+lgb[l]])
            end
            temp = [reduce_mirror(gbasis[l][i][1], N), reduce_mirror(gbasis[l][j][1], N)]
            Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; temp]))
            if i == j
                @inbounds add_to_expression!(cons[Locb], gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], 2*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            end
        end
    end
    mons = generate_mons(N, lso)
    mons = filter_mons(mons, tsupp, H, N, identify_zeros=false)
    for mon in mons
        fr = @variable(model)
        for (i, item) in enumerate(H.supp)
            bi,c = reduce!([item; mon], N, identify_zeros=false)
            if imag(c) != 0
                Locb = bfind(tsupp, [bi])
                add_to_expression!(cons[Locb], H.coe[i]*imag(c), fr)
            end
        end
    end
    # posepsd7!(model, cons, tsupp, N, identify_zeros=false)
    @variable(model, λ)
    cons[1] += λ
    for i = 1:N
        μ = @variable(model, lower_bound = 0)
        add_to_expression!(cons[1], μ)
        for j = 1:3
            Locb = bfind(tsupp, [[3(i-1)+j], [3(i-1)+j]])
            add_to_expression!(cons[Locb], -μ)
        end
    end
    @objective(model, Max, λ)
    @constraint(model, cons .== 0)
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
    flag = status == MathOptInterface.OPTIMAL ? 1 : 0
    return flag
end

function certify_Heisenberg_kagome_gap(N::Int, H::ncpoly, triples, edges, inner_triples, inner_edges, gamma, d::Int; lso=4, QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    basis = [get_kagome_basis(N, triples, edges, d, label=i) for i = 1:2]
    lb = length.(basis)
    gbasis = [get_kagome_bulkbasis(N, inner_triples, inner_edges, d-1, label=i) for i = 1:2]
    lgb = length.(gbasis)
    bs = [lb; lgb]
    if QUIET == false
        println("The block sizes are $bs.")
    end
    tsupp = Vector{Vector{Int}}[]
    for i = 1:length(basis), j = 1:lb[i], k = j:lb[i]
        @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]], N, model="kagome")
        if c != 0
            if isempty(bi)
                push!(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
            else
                push!(tsupp, sort([basis[i][j][2]; basis[i][k][2]; [bi]]))
            end
        end
    end
    for l = 1:length(gbasis), i = 1:lgb[l], j = i:lgb[l]
        bis = PSDstate_entry(gbasis[l][i][1], gbasis[l][j][1], H, N, model="kagome")[1]
        for bi in bis
            if isempty(bi)
                push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
            else
                push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
            end
        end
        if !isz(gbasis[l][i][1], model="kagome") && !isz(gbasis[l][j][1], model="kagome")
            temp = [reduce_perm(gbasis[l][i][1]), reduce_perm(gbasis[l][j][1])]
            push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; temp]))
        end
    end
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3, w = 0:3
        ind = [i,j,k,l,s,t,u,v,w]
        if all(x->iseven(sum(ind .== x)), 1:3)
            inx = ind .!= 0
            push!(tsupp, [reduce_perm(3*(Vector(1:9)[inx] .- 1) + ind[inx])])
        end
    end
    sort!(tsupp)
    unique!(tsupp)
    if QUIET == false
        println("There are $(length(tsupp)) affine constraints.")
        println("Assembling the SDP...")
    end
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:length(tsupp)]
    pos = Vector{Symmetric{VariableRef}}(undef, length(basis))
    for i = 1:length(basis)
        pos[i] = @variable(model, [1:lb[i], 1:lb[i]], PSD)
        for j = 1:lb[i], k = j:lb[i]
            @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]], N, realify=true, model="kagome")
            if c != 0
                if isempty(bi)
                    Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
                else
                    Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]; [bi]]))
                end
                if j == k
                    @inbounds add_to_expression!(cons[Locb], c, pos[i][j,k])
                else
                    @inbounds add_to_expression!(cons[Locb], 2c, pos[i][j,k])
                end
            end
        end
    end
    gpos = Vector{Symmetric{VariableRef}}(undef, length(gbasis))
    for l = 1:length(gbasis)
        gpos[l] = @variable(model, [1:lgb[l], 1:lgb[l]], PSD)
        for i = 1:lgb[l], j = i:lgb[l]
            bis,coes = PSDstate_entry(gbasis[l][i][1], gbasis[l][j][1], H, N, realify=true, model="kagome")
            for (k, bi) in enumerate(bis)
                if isempty(bi)
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
                else
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
                end
                if i == j
                    @inbounds add_to_expression!(cons[Locb], coes[k], gpos[l][i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], 2*coes[k], gpos[l][i, j])
                end
            end
            @inbounds bi,c = reduce!([gbasis[l][i][1]; gbasis[l][j][1]], N, realify=true, model="kagome")
            if c != 0
                if isempty(bi)
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
                else
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
                end
                if i == j
                    @inbounds add_to_expression!(cons[Locb], -c*gamma, gpos[l][i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], -2c*gamma, gpos[l][i, j])
                end
            end
            if !isz(gbasis[l][i][1], model="kagome") && !isz(gbasis[l][j][1], model="kagome")
                temp = [reduce_perm(gbasis[l][i][1]), reduce_perm(gbasis[l][j][1])]
                Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; temp]))
                if i == j
                    @inbounds add_to_expression!(cons[Locb], gamma, gpos[l][i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], 2*gamma, gpos[l][i, j])
                end
            end
        end
    end
    if N > 5
        if N == 13
            iN = 5
        elseif N == 27
            iN = 13
        elseif N == 45
            iN = 27
        else
            @error "Wrong number of sites!"
        end
        mons = generate_mons_kagome(iN, lso)
        mons = filter_mons(mons, tsupp, H, N, model="kagome")
        for mon in mons
            fr = @variable(model)
            for (i, item) in enumerate(H.supp)
                bi,c = reduce!([item; mon], N, model="kagome")
                if imag(c) != 0
                    Locb = bfind(tsupp, [bi])
                    add_to_expression!(cons[Locb], H.coe[i]*imag(c), fr)
                end
            end
        end
    end
    posepsd9!(model, cons, tsupp)
    @variable(model, λ)
    cons[1] += λ
    @objective(model, Max, λ)
    @constraint(model, cons .== 0)
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
    flag = status == MathOptInterface.OPTIMAL ? 1 : 0
    return flag
end

function certify_Heisenberg_kagome_gap_nosignsymmetry(N::Int, H::ncpoly, triples, edges, inner_triples, inner_edges, gamma, d::Int; obj=nothing, lso=4, QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    # basis = [get_kagome_basis(N, triples, edges, d, label=i) for i = 1:4]
    basis = [vcat([get_kagome_basis(N, triples, edges, d, label=i) for i = 1:4]...)]
    lb = length.(basis)
    # gbasis = [get_kagome_bulkbasis(N, inner_triples, inner_edges, d-1, label=i) for i = 2:4]
    gbasis = [vcat([get_kagome_bulkbasis(N, inner_triples, inner_edges, d-1, label=i) for i = 1:4]...)]
    lgb = length.(gbasis)
    bs = [lb; lgb]
    if QUIET == false
        println("The block sizes are $bs.")
    end
    tsupp = Vector{Vector{Int}}[]
    for i = 1:length(basis), j = 1:lb[i], k = j:lb[i]
        @inbounds bi = reduce!([basis[i][j][1]; basis[i][k][1]], N, identify_zeros=false, symmetry=false)[1]
        if isempty(bi)
            push!(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
        else
            push!(tsupp, sort([basis[i][j][2]; basis[i][k][2]; [bi]]))
        end
        temp = [basis[i][j][1], basis[i][k][1]]
        push!(tsupp, sort([basis[i][j][2]; basis[i][k][2]; temp]))
    end
    for l = 1:length(gbasis), i = 1:lgb[l], j = i:lgb[l]
        bis = PSDstate_entry(gbasis[l][i][1], gbasis[l][j][1], H, N, identify_zeros=false, symmetry=false)[1]
        for bi in bis
            if isempty(bi)
                push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
            else
                push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
            end
        end
        temp = [gbasis[l][i][1], gbasis[l][j][1]]
        push!(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; temp]))
    end
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3
        ind = [i,j,k,l,s,t,u]
        inx = ind .!= 0
        push!(tsupp, [reduce_perm(3*(Vector(1:7)[inx] .- 1) + ind[inx])])
    end
    for i = 1:N, j = 1:3
        push!(tsupp, [[3(i-1)+j], [3(i-1)+j]])
    end
    sort!(tsupp)
    unique!(tsupp)
    if QUIET == false
        println("There are $(length(tsupp)) affine constraints.")
        println("Assembling the SDP...")
    end
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:length(tsupp)]
    pos = Vector{Symmetric{VariableRef}}(undef, length(basis))
    for i = 1:length(basis)
        pos[i] = @variable(model, [1:2*lb[i], 1:2*lb[i]], PSD)
        for j = 1:lb[i], k = j:lb[i]
            @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]], N, identify_zeros=false, symmetry=false)
            if isempty(bi)
                Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
            else
                Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]; [bi]]))
            end
            if j == k
                @inbounds add_to_expression!(cons[Locb], c, pos[i][j,k]+pos[i][j+lb[i],k+lb[i]])
            elseif isreal(c)
                @inbounds add_to_expression!(cons[Locb], 2c, pos[i][j,k]+pos[i][j+lb[i],k+lb[i]])
            else
                @inbounds add_to_expression!(cons[Locb], -2*imag(c), pos[i][j,k+lb[i]]-pos[i][k,j+lb[i]])
            end
        end
    end
    npos = Vector{Symmetric{VariableRef}}(undef, length(basis))
    for l = 1:length(basis)
        npos[l] = @variable(model, [1:2*lb[l], 1:2*lb[l]], PSD)
        for i = 1:lb[l], j = i:lb[l]
            @inbounds bi,c = reduce!([basis[l][i][1]; basis[l][j][1]], N, identify_zeros=true, symmetry=true)
            if isempty(bi)
                Locb = bfind(tsupp, sort([basis[l][i][2]; basis[l][j][2]]))
            else
                Locb = bfind(tsupp, sort([basis[l][i][2]; basis[l][j][2]; [bi]]))
            end
            if i == j
                @inbounds add_to_expression!(cons[Locb], c, npos[l][i, j]+npos[l][i+lb[l], j+lb[l]])
            elseif isreal(c)
                @inbounds add_to_expression!(cons[Locb], 2c, npos[l][i, j]+npos[l][i+lb[l], j+lb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], -2*imag(c), npos[l][i,j+lb[l]]-npos[l][j,i+lb[l]])
            end
            temp = [basis[l][i][1], basis[l][j][1]]
            Locb = bfind(tsupp, sort([basis[l][i][2]; basis[l][j][2]; temp]))
            if i == j
                @inbounds add_to_expression!(cons[Locb], -1, npos[l][i, j]+npos[l][i+lb[l], j+lb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], -2, npos[l][i, j]+npos[l][i+lb[l], j+lb[l]])
            end
        end
    end
    gpos = Vector{Symmetric{VariableRef}}(undef, length(gbasis))
    for l = 1:length(gbasis)
        gpos[l] = @variable(model, [1:2*lgb[l], 1:2*lgb[l]], PSD)
        for i = 1:lgb[l], j = i:lgb[l]
            bis,coes = PSDstate_entry(gbasis[l][i][1], gbasis[l][j][1], H, N, identify_zeros=false, symmetry=false)
            for (k, bi) in enumerate(bis)
                if isempty(bi)
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
                else
                    Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
                end
                if i == j
                    @inbounds add_to_expression!(cons[Locb], coes[k], gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
                elseif isreal(coes[k])
                    @inbounds add_to_expression!(cons[Locb], 2*coes[k], gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
                else
                    @inbounds add_to_expression!(cons[Locb], -2*imag(coes[k]), gpos[l][i, j+lgb[l]]-gpos[l][j, i+lgb[l]])
                end
            end
            @inbounds bi,c = reduce!([gbasis[l][i][1]; gbasis[l][j][1]], N, identify_zeros=false, symmetry=false)
            if isempty(bi)
                Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]]))
            else
                Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; [bi]]))
            end
            if i == j
                @inbounds add_to_expression!(cons[Locb], -c*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            elseif isreal(c)
                @inbounds add_to_expression!(cons[Locb], -2c*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], 2*imag(c)*gamma, gpos[l][i, j+lgb[l]]-gpos[l][j, i+lgb[l]])
            end
            temp = [gbasis[l][i][1], gbasis[l][j][1]]
            Locb = bfind(tsupp, sort([gbasis[l][i][2]; gbasis[l][j][2]; temp]))
            if i == j
                @inbounds add_to_expression!(cons[Locb], gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], 2*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            end
        end
    end
    if N > 5
        if N == 13
            iN = 5
        elseif N == 27
            iN = 13
        elseif N == 45
            iN = 27
        else
            @error "Wrong number of sites!"
        end
        mons = generate_mons_kagome(iN, lso)
        mons = filter_mons(mons, tsupp, H, N, model="kagome", identify_zeros=false, symmetry=false)
        for mon in mons
            fr = @variable(model)
            for (i, item) in enumerate(H.supp)
                bi,c = reduce!([item; mon], N, model="kagome", identify_zeros=false, symmetry=false)
                if imag(c) != 0
                    Locb = bfind(tsupp, [bi])
                    add_to_expression!(cons[Locb], H.coe[i]*imag(c), fr)
                end
            end
        end
    end
    posepsd7!(model, cons, tsupp, N, lattice="kagome", identify_zeros=true, symmetry=true)
    @variable(model, λ)
    cons[1] -= λ
    if obj !== nothing
        for i = 1:length(obj.supp)
            Locb = bfind(tsupp, [obj.supp[i]])
            cons[Locb] += obj.coe[i]
        end
    end
    for i = 1:N
        μ = @variable(model, lower_bound = 0)
        add_to_expression!(cons[1], μ)
        for j = 1:3
            Locb = bfind(tsupp, [[3(i-1)+j], [3(i-1)+j]])
            add_to_expression!(cons[Locb], -μ)
        end
    end
    @objective(model, Min, λ)
    @constraint(model, cons .== 0)
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
    flag = status == MathOptInterface.OPTIMAL ? 1 : 0
    return flag,value(λ)
end
