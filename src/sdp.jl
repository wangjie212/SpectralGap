function certify_gap(N::Int, H::ncpoly, gamma, d::Int; QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    basis = [get_basis(N, d, label=i) for i in [1,2]]
    lb = length.(basis)
    gbasis = [get_wbasis(N, d-1, label=i) for i in [1,2]] 
    lgb = length.(gbasis)
    if d > 2
        tbasis = [tuple(a, b) for a in get_wbasis(N, 1, label=1), b in get_sbasis(N)]
        ltb = length(tbasis)
    end
    if QUIET == false
        println("The block sizes are $lb.")
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
            if !isz(gbasis[l][i]) && !isz(gbasis[l][j])
                Locb = bfind(tsupp, sort([reduce4(gbasis[l][i], N), reduce4(gbasis[l][j], N)]))
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
            if !isz(tbasis[i][1]) && !isz(tbasis[j][1])
                Locb = bfind(tsupp, sort([tbasis[i][2]; tbasis[j][2]; [reduce4(tbasis[i][1], N)]; [reduce4(tbasis[j][1], N)]]))
                if i == j
                    @inbounds add_to_expression!(cons[Locb], gamma, tpos[i, j])
                else
                    @inbounds add_to_expression!(cons[Locb], 2*gamma, tpos[i, j])
                end
            end
        end
    end
    @variable(model, lower)
    cons[1] += lower
    # @constraint(model, cons .== 0)
    @constraint(model, con, cons==zeros(length(cons)))
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
    dual_var = -dual(con)
    mmat = Vector{Matrix{ComplexF16}}(undef, length(basis))
    for i = 1:length(basis)
        mmat[i] = zeros(ComplexF16, lb[i], lb[i])
        for j = 1:lb[i], k = 1:lb[i]
            @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]], N)
            if c != 0
                if isempty(bi)
                    if isempty(basis[i][j][2]) && isempty(basis[i][k][2])
                        Locb = 1
                    else
                        Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]]))
                    end
                else
                    Locb = bfind(tsupp, sort([basis[i][j][2]; basis[i][k][2]; [bi]]))
                end
                mmat[i][j,k] = c*dual_var[Locb]
            end
        end
    end
    return status,basis,mmat,tsupp
end
