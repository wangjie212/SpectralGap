function certify_gap(N::Int, H::ncpoly, gamma, d::Int; QUIET=false)
    println("********************************** SpectralGap **********************************")
    println("SpectralGap is launching...")
    basis = [get_basis(N, d, label=i) for i = 1:3]
    lb = length.(basis)
    gbasis = get_wbasis(N, d-1)
    lgb = length(gbasis)
    if QUIET == false
        println("The block sizes are $lb.")
    end
    tsupp = [[Int[]]]
    for i = 1:3, j = 1:lb[i], k = j:lb[i]
        @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]])
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
    tsupp = tsupp[2:end]
    if QUIET == false
        println("There are $(length(tsupp)) affine constraints.")
        println("Assembling the SDP...")
    end
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)
    cons = [AffExpr(0) for i=1:length(tsupp)]
    pos = Vector{Symmetric{VariableRef}}(undef, 3)
    for i = 1:3
        pos[i] = @variable(model, [1:lb[i], 1:lb[i]], PSD)
        for j = 1:lb[i], k = j:lb[i]
            @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]], realify=true)
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
                if j == k
                    @inbounds add_to_expression!(cons[Locb], c, pos[i][j,k])
                else
                    @inbounds add_to_expression!(cons[Locb], 2c, pos[i][j,k])
                end
            end
        end
    end
    gpos = @variable(model, [1:2*lgb, 1:2*lgb], PSD)
    for i = 1:lgb, j = i:lgb
        pp1 = gpos[i, j] + gpos[i+lgb, j+lgb]
        pp2 = gpos[i+lgb, j] - gpos[j+lgb, i]
        for k = 1:length(H.supp)
            @inbounds bi,c = reduce!([gbasis[i]; H.supp[k]; gbasis[j]])
            if c != 0
                Locb = bfind(tsupp, [bi])
                println(bi)
                if i == j
                    @inbounds add_to_expression!(cons[Locb], real(c)*H.coe[k], pp1)
                else
                    @inbounds add_to_expression!(cons[Locb], 2*real(c)*H.coe[k], pp1)
                    if imag(c) != 0
                        @inbounds add_to_expression!(cons[Locb], -2*imag(c)*H.coe[k], pp2)
                    end
                end
            end
            @inbounds bi,c = reduce!([gbasis[i]; gbasis[j]; H.supp[k]])
            if c != 0
                Locb = bfind(tsupp, [bi])
                if i == j
                    @inbounds add_to_expression!(cons[Locb], -real(c)*H.coe[k], pp1)
                else
                    @inbounds add_to_expression!(cons[Locb], -2*real(c)*H.coe[k], pp1)
                    if imag(c) != 0
                        @inbounds add_to_expression!(cons[Locb], 2*imag(c)*H.coe[k], pp2)
                   end
                end
            end
        end
        @inbounds bi,c = reduce!([gbasis[i]; gbasis[j]])
        if c != 0
            Locb = bfind(tsupp, [bi])
            if i == j
                @inbounds add_to_expression!(cons[Locb], -gamma*real(c), pp1)
            else
                @inbounds add_to_expression!(cons[Locb], -2*gamma*real(c), pp1)
                if imag(c) != 0
                    @inbounds add_to_expression!(cons[Locb], 2*gamma*imag(c), pp2)
                end
            end
        end
        if !isz(gbasis[i]) && !isz(gbasis[j])
            Locb = bfind(tsupp, sort([gbasis[i], gbasis[j]]))
            if i == j
                @inbounds add_to_expression!(cons[Locb], gamma, pp1)
            else
                @inbounds add_to_expression!(cons[Locb], 2*gamma, pp1)
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
    mmat = Vector{Matrix{ComplexF16}}(undef, 3)
    for i = 1:3
        mmat[i] = zeros(ComplexF16, lb[i], lb[i])
        for j = 1:lb[i], k = 1:lb[i]
            @inbounds bi,c = reduce!([basis[i][j][1]; basis[i][k][1]])
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
    return status,basis,mmat
end
