function certify_Ising_gap(N::Int, H::ncpoly, gamma, d::Int; QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    basis = [get_basis(N, d, label=i) for i in [1]]
    lb = length.(basis)
    gbasis = [get_wbasis(N, d-1, label=i) for i in [1,2]] 
    lgb = length.(gbasis)
    bs = [lb; lgb]
    if d > 2
        tbasis = [tuple(a, b) for a in get_wbasis(N, 1, label=1), b in get_sbasis(N, 1)]
        ltb = length(tbasis)
        bs = [bs; ltb]
    end
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
        # for k = 1:length(H.supp)
        #     if !isreal(reduce2!(reduce3!(reduce1!([H.supp[k]; gbasis[l][j]])))[2])
        #         @inbounds bi,c = reduce!([gbasis[l][i]; H.supp[k]; gbasis[l][j]], N)
        #         if c != 0 && !isempty(bi)
        #             push!(tsupp, [bi])
        #         end
        #     end
        # end
        # @inbounds bi,c = reduce!([gbasis[l][i]; gbasis[l][j]], N)
        # if c != 0 && !isempty(bi)
        #     push!(tsupp, [bi])
        # end
        if !isz(gbasis[l][i], model="Ising") && !isz(gbasis[l][j], model="Ising")
            push!(tsupp, sort([reduce_mirror(gbasis[l][i], N), reduce_mirror(gbasis[l][j], N)]))
        end
    end
    if d > 2
        for i = 1:ltb, j = i:ltb
            for k = 1:length(H.supp)
                if !isreal(reduce2!(reduce3!(reduce1!([H.supp[k]; tbasis[j][1]])))[2])
                    @inbounds bi,c = reduce!([tbasis[i][1]; H.supp[k]; tbasis[j][1]], N)
                    if c != 0
                        if isempty(bi)
                            push!(tsupp, sort([tbasis[i][2]; tbasis[j][2]]))
                        else
                            push!(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [bi]]))
                        end
                    end
                end
            end
            @inbounds bi,c = reduce!([tbasis[i][1]; tbasis[j][1]], N)
            if c != 0
                if isempty(bi)
                    push!(tsupp, sort([tbasis[i][2]; tbasis[j][2]]))
                else
                    push!(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [bi]]))
                end
            end
            if !isz(tbasis[i][1], model="Ising") && !isz(tbasis[j][1], model="Ising")
                push!(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [reduce_mirror(tbasis[i][1], N)]; [reduce_mirror(tbasis[j][1], N)]]))
            end
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
            for k = 1:length(H.supp)
                if !isreal(reduce2!(reduce3!(reduce1!([H.supp[k]; gbasis[l][j]])))[2])
                    @inbounds bi,c = reduce!([gbasis[l][i]; H.supp[k]; gbasis[l][j]], N, realify=true)
                    if c != 0
                        Locb = isempty(bi) ? 1 : bfind(tsupp, [bi])
                        if i == j
                            @inbounds add_to_expression!(cons[Locb], 2c*H.coe[k], gpos[l][i, j])
                        else
                            @inbounds add_to_expression!(cons[Locb], 4c*H.coe[k], gpos[l][i, j])
                        end
                    end
                end
            end
            @inbounds bi,c = reduce!([gbasis[l][i]; gbasis[l][j]], N, realify=true)
            if c != 0
                Locb = isempty(bi) ? 1 : bfind(tsupp, [bi])
                if i == j
                    @inbounds add_to_expression!(cons[Locb], -c*gamma, gpos[l][i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], -2c*gamma, gpos[l][i, j])
                end
            end
            if !isz(gbasis[l][i], model="Ising") && !isz(gbasis[l][j], model="Ising")
                Locb = bfind(tsupp, sort([reduce_mirror(gbasis[l][i], N), reduce_mirror(gbasis[l][j], N)]))
                if i == j
                    @inbounds add_to_expression!(cons[Locb], gamma, gpos[l][i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], 2*gamma, gpos[l][i, j])
                end
            end
        end
    end
    if d > 2
        tpos = @variable(model, [1:ltb, 1:ltb], PSD)
        for i = 1:ltb, j = i:ltb
            for k = 1:length(H.supp)
                if !isreal(reduce2!(reduce3!(reduce1!([H.supp[k]; tbasis[j][1]])))[2])
                    @inbounds bi,c = reduce!([tbasis[i][1]; H.supp[k]; tbasis[j][1]], N, realify=true)
                    if c != 0
                        if isempty(bi)
                            Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]]))
                        else
                            Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [bi]]))
                        end
                        if i == j
                            @inbounds add_to_expression!(cons[Locb], 2c*H.coe[k], tpos[i, j])
                        else
                            @inbounds add_to_expression!(cons[Locb], 4c*H.coe[k], tpos[i, j])
                        end
                    end
                end
            end
            @inbounds bi,c = reduce!([tbasis[i][1]; tbasis[j][1]], N, realify=true)
            if c != 0
                if isempty(bi)
                    Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]]))
                else
                    Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [bi]]))
                end
                if i == j
                    @inbounds add_to_expression!(cons[Locb], -c*gamma, tpos[i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], -2c*gamma, tpos[i, j])
                end
            end
            if !isz(tbasis[i][1], model="Ising") && !isz(tbasis[j][1], model="Ising")
                Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [reduce_mirror(tbasis[i][1], N)]; [reduce_mirror(tbasis[j][1], N)]]))
                if i == j
                    @inbounds add_to_expression!(cons[Locb], gamma, tpos[i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], 2*gamma, tpos[i, j])
                end
            end
        end
    end
    @variable(model, λ)
    cons[1] += λ
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
    #                 if isempty(basis[i][j][2]) && isempty(basis[i][k][2])
    #                     Locb = 1
    #                 else
    #                     Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
    #                 end
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

function certify_Ising_gap_nosignsymmetry(N::Int, H::ncpoly, gamma, d::Int; QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    basis = [get_basis(N, d, label=1)]
    lb = length.(basis)
    gbasis = [[get_wbasis(N, d-1, label=1); get_wbasis(N, d-1, label=2)]]
    lgb = length.(gbasis)
    bs = [lb; lgb]
    if d > 2
        tbasis = [tuple(a, b) for a in [get_wbasis(N, 1, label=1); get_wbasis(N, 1, label=2)], b in get_sbasis(N, 1)]
        ltb = length(tbasis)
        bs = [bs; ltb]
    end
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
    end
    for l = 1:length(gbasis), i = 1:lgb[l], j = i:lgb[l]
        for k = 1:length(H.supp)
            if !isreal(reduce2!(reduce3!(reduce1!([H.supp[k]; gbasis[l][j]])))[2])
                @inbounds bi = reduce!([gbasis[l][i]; H.supp[k]; gbasis[l][j]], N, identify_zeros=false)[1]
                push!(tsupp, [bi])
            end
        end
        # @inbounds bi = reduce!([gbasis[l][i]; gbasis[l][j]], N, identify_zeros=false)[1]
        # if !isempty(bi)
        #     push!(tsupp, [bi])
        # end
        push!(tsupp, sort([reduce_mirror(gbasis[l][i], N), reduce_mirror(gbasis[l][j], N)]))
    end
    if d > 2
        for i = 1:ltb, j = i:ltb
            for k = 1:length(H.supp)
                if !isreal(reduce2!(reduce3!(reduce1!([H.supp[k]; tbasis[j][1]])))[2])
                    @inbounds bi = reduce!([tbasis[i][1]; H.supp[k]; tbasis[j][1]], N, identify_zeros=false)[1]
                    if isempty(bi)
                        push!(tsupp, sort([tbasis[i][2]; tbasis[j][2]]))
                    else
                        push!(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [bi]]))
                    end
                end
            end
            @inbounds bi = reduce!([tbasis[i][1]; tbasis[j][1]], N, identify_zeros=false)[1]
            if isempty(bi)
                push!(tsupp, sort([tbasis[i][2]; tbasis[j][2]]))
            else
                push!(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [bi]]))
            end
            push!(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [reduce_mirror(tbasis[i][1], N)]; [reduce_mirror(tbasis[j][1], N)]]))
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
    gpos = Vector{Symmetric{VariableRef}}(undef, length(gbasis))
    for l = 1:length(gbasis)
        gpos[l] = @variable(model, [1:2*lgb[l], 1:2*lgb[l]], PSD)
        for i = 1:lgb[l], j = i:lgb[l]
            for k = 1:length(H.supp)
                if !isreal(reduce2!(reduce3!(reduce1!([H.supp[k]; gbasis[l][j]])))[2])
                    @inbounds bi,c = reduce!([gbasis[l][i]; H.supp[k]; gbasis[l][j]], N, identify_zeros=false)
                    Locb = isempty(bi) ? 1 : bfind(tsupp, [bi])
                    if i == j
                        @inbounds add_to_expression!(cons[Locb], 2c*H.coe[k], gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
                    elseif isreal(c)
                        @inbounds add_to_expression!(cons[Locb], 4c*H.coe[k], gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
                    else
                        @inbounds add_to_expression!(cons[Locb], -4*imag(c)*H.coe[k], gpos[l][i, j+lgb[l]]-gpos[l][j, i+lgb[l]])
                    end
                end
            end
            @inbounds bi,c = reduce!([gbasis[l][i]; gbasis[l][j]], N, identify_zeros=false)
            Locb = isempty(bi) ? 1 : bfind(tsupp, [bi])
            if i == j
                @inbounds add_to_expression!(cons[Locb], -c*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            elseif isreal(c)
                @inbounds add_to_expression!(cons[Locb], -2c*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], 2*imag(c)*gamma, gpos[l][i, j+lgb[l]]-gpos[l][j, i+lgb[l]])
            end
            Locb = bfind(tsupp, sort([reduce_mirror(gbasis[l][i], N), reduce_mirror(gbasis[l][j], N)]))
            if i == j
                @inbounds add_to_expression!(cons[Locb], gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], 2*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            end
        end
    end
    if d > 2
        tpos = @variable(model, [1:2*ltb, 1:2*ltb], PSD)
        for i = 1:ltb, j = i:ltb
            for k = 1:length(H.supp)
                if !isreal(reduce2!(reduce3!(reduce1!([H.supp[k]; tbasis[j][1]])))[2])
                    @inbounds bi,c = reduce!([tbasis[i][1]; H.supp[k]; tbasis[j][1]], N, identify_zeros=false)
                    if isempty(bi)
                        Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]]))
                    else
                        Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [bi]]))
                    end
                    if i == j
                        @inbounds add_to_expression!(cons[Locb], 2c*H.coe[k], tpos[i, j]+tpos[i+ltb, j+ltb])
                    elseif isreal(c)
                        @inbounds add_to_expression!(cons[Locb], 4c*H.coe[k], tpos[i, j]+tpos[i+ltb, j+ltb])
                    else
                        @inbounds add_to_expression!(cons[Locb], -4*imag(c)*H.coe[k], tpos[i, j+ltb]-tpos[j, i+ltb])
                    end
                end
            end
            @inbounds bi,c = reduce!([tbasis[i][1]; tbasis[j][1]], N, identify_zeros=false)
            if isempty(bi)
                Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]]))
            else
                Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [bi]]))
            end
            if i == j
                @inbounds add_to_expression!(cons[Locb], -c*gamma, tpos[i, j]+tpos[i+ltb, j+ltb])
            elseif isreal(c)
                @inbounds add_to_expression!(cons[Locb], -2c*gamma, tpos[i, j]+tpos[i+ltb, j+ltb])
            else
                @inbounds add_to_expression!(cons[Locb], 2*imag(c)*gamma, tpos[i, j+ltb]-tpos[j, i+ltb])
            end
            Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [reduce_mirror(tbasis[i][1], N)]; [reduce_mirror(tbasis[j][1], N)]]))
            if i == j
                @inbounds add_to_expression!(cons[Locb], gamma, tpos[i, j]+tpos[i+ltb, j+ltb])
            else
                @inbounds add_to_expression!(cons[Locb], 2*gamma, tpos[i, j]+tpos[i+ltb, j+ltb])
            end
        end
    end
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

function certify_Heisenberg_kagome_gap(N::Int, H::ncpoly, triples, edges, inner_triples, inner_edges, gamma, d::Int; QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    basis = [get_kagome_basis(N, triples, edges, d, label=i) for i = 1:2]
    lb = length.(basis)
    gbasis = [get_kagome_wbasis(N, inner_triples, inner_edges, d-1, label=i) for i = 1:2]
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
        if !isz(gbasis[l][i], model="kagome") && !isz(gbasis[l][j], model="kagome")
            push!(tsupp, sort([reduce_perm(gbasis[l][i]), reduce_perm(gbasis[l][j])]))
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
            for k = 1:length(H.supp)
                if !isreal(reduce2!(reduce3!(reduce1!([H.supp[k]; gbasis[l][j]])))[2])
                    @inbounds bi,c = reduce!([gbasis[l][i]; H.supp[k]; gbasis[l][j]], N, realify=true, model="kagome")
                    if c != 0
                        Locb = isempty(bi) ? 1 : bfind(tsupp, [bi])
                        if i == j
                            @inbounds add_to_expression!(cons[Locb], c*H.coe[k]/2, gpos[l][i, j])
                        else
                            @inbounds add_to_expression!(cons[Locb], c*H.coe[k], gpos[l][i, j])
                        end
                    end
                end
            end
            @inbounds bi,c = reduce!([gbasis[l][i]; gbasis[l][j]], N, realify=true, model="kagome")
            if c != 0
                Locb = isempty(bi) ? 1 : bfind(tsupp, [bi])
                if i == j
                    @inbounds add_to_expression!(cons[Locb], -c*gamma, gpos[l][i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], -2c*gamma, gpos[l][i, j])
                end
            end
            if !isz(gbasis[l][i], model="kagome") && !isz(gbasis[l][j], model="kagome")
                Locb = bfind(tsupp, sort([reduce_perm(gbasis[l][i]), reduce_perm(gbasis[l][j])]))
                if i == j
                    @inbounds add_to_expression!(cons[Locb], gamma, gpos[l][i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], 2*gamma, gpos[l][i, j])
                end
            end
        end
    end
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

function certify_Heisenberg_kagome_gap_nosignsymmetry(N::Int, H::ncpoly, triples, edges, gamma, d::Int; obj=nothing, QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    basis = [get_kagome_basis(N, triples, edges, d, label=i) for i = 1:4]
    lb = length.(basis)
    gbasis = [get_kagome_wbasis(N, [], [], d-1, label=i) for i = 2:4]
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
    end
    # for l = 1:length(gbasis), i = 1:lgb[l], j = i:lgb[l]
    #     push!(tsupp, sort([gbasis[l][i], gbasis[l][j]]))
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
    gpos = Vector{Symmetric{VariableRef}}(undef, length(gbasis))
    for l = 1:length(gbasis)
        gpos[l] = @variable(model, [1:2*lgb[l], 1:2*lgb[l]], PSD)
        for i = 1:lgb[l], j = i:lgb[l]
            for k = 1:length(H.supp)
                if !isreal(reduce2!(reduce3!(reduce1!([H.supp[k]; gbasis[l][j]])))[2])
                    @inbounds bi,c = reduce!([gbasis[l][i]; H.supp[k]; gbasis[l][j]], N, identify_zeros=false, symmetry=false)
                    Locb = isempty(bi) ? 1 : bfind(tsupp, [bi])
                    if i == j
                        @inbounds add_to_expression!(cons[Locb], c*H.coe[k]/2, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
                    elseif isreal(c)
                        @inbounds add_to_expression!(cons[Locb], c*H.coe[k], gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
                    else
                        @inbounds add_to_expression!(cons[Locb], -imag(c)*H.coe[k], gpos[l][i, j+lgb[l]]-gpos[l][j, i+lgb[l]])
                    end
                end
            end
            @inbounds bi,c = reduce!([gbasis[l][i]; gbasis[l][j]], N, identify_zeros=false, symmetry=false)
            Locb = isempty(bi) ? 1 : bfind(tsupp, [bi])
            if i == j
                @inbounds add_to_expression!(cons[Locb], -c*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            elseif isreal(c)
                @inbounds add_to_expression!(cons[Locb], -2c*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            else
                @inbounds add_to_expression!(cons[Locb], 2imag(c)*gamma, gpos[l][i, j+lgb[l]]-gpos[l][j, i+lgb[l]])
            end
            # Locb = bfind(tsupp, sort([gbasis[l][i], gbasis[l][j]]))
            # if i == j
            #     @inbounds add_to_expression!(cons[Locb], gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            # else
            #     @inbounds add_to_expression!(cons[Locb], 2*gamma, gpos[l][i, j]+gpos[l][i+lgb[l], j+lgb[l]])
            # end
        end
    end
    @variable(model, λ)
    cons[1] -= λ
    if obj !== nothing
        for i = 1:length(obj.supp)
            Locb = bfind(tsupp, [obj.supp[i]])
            cons[Locb] += obj.coe[i]
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
