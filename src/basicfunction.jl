mutable struct ncpoly
    supp::Vector{Vector{Int}}
    coe::Vector{Float64}
end

# reduction to the normal form
function reduce1!(a::Vector{Int})
    la = length(a)
    flag = 1
    while flag == 1
        ind = findfirst(x->ceil(Int, a[x]/3) > ceil(Int, a[x+1]/3), 1:la-1)
        if ind !== nothing
            a[ind],a[ind+1] = a[ind+1],a[ind]
            flag = 1
        else
            flag = 0
        end
    end
    return a
end

# reduction to the normal form
function reduce2!(a::Vector{Int}; realify=false)
    la = length(a)
    flag = 1
    coef = 1
    while flag == 1
        ind = findfirst(x -> a[x] != a[x+1] && ceil(Int, a[x]/3) == ceil(Int, a[x+1]/3), 1:la-1)
        if ind !== nothing
            s = mod.(a[ind:ind+1], 3)
            if s == [1, 2]
                a[ind] += Int(2)
                coef *= im
            elseif s == [0, 2]
                a[ind] -= Int(2)
                coef *= -im
            elseif s == [2, 1] || s == [1, 0]
                a[ind] += Int(1)
                coef *= -im
            else
                a[ind] -= Int(1)
                coef *= im  
            end
            deleteat!(a, ind+1)
            la -= 1
            flag = 1
        else
            flag = 0
        end
    end
    if realify == true
        coef = isreal(coef) ? real(coef) : imag(coef)
    end
    return a,coef
end

# reduction to the normal form
function reduce3!(a::Vector{Int})
    i = 1
    while i < length(a)
        if a[i] == a[i+1]
            deleteat!(a, i)
            deleteat!(a, i)
        else
            i += 1
        end
    end
    return a
end

# identify zeros by sign symmetry
function isz(a::Vector{Int}; model="Ising")
    if model == "Ising"
        return any(i->isodd(count(isequal(i), mod.(a,3))), [0,2]) # for Ising model
    else
        return any(i->isodd(count(isequal(i), mod.(a,3))), [0,1,2]) # for Kagome lattice
    end
end

# reduction w.r.t permutation symmetry
function reduce_perm(a::Vector{Int})
    ra = smod.(a, 3)
    # sym = [[1;2;3], [1;3;2], [2;1;3], [2;3;1], [3;1;2], [3;2;1]]
    sym = [[1;2;3], [2;3;1], [3;1;2]]
    return findmin([Int.(3*(ceil.(Int, a./3).-1) .+ s[ra]) for s in sym])[1]
end

# reduction w.r.t mirror symmetry
function reduce_mirror(a::Vector{Int}, N)
    ra = reverse(a)
    ma = 3*(N .- ceil.(Int, ra/3)) + smod.(ra, 3)
    return min(a, ma)
end

function smod(i, s)
    r = mod(i, s)
    return r == 0 ? s : r
end

# implement all reductions
function reduce!(a::Vector{Int}, N; model="Ising", realify=false, identify_zeros=true, symmetry=true)
    reduce1!(a)
    reduce3!(a)
    a,coef = reduce2!(a, realify=realify)
    reduce3!(a)
    if !isempty(a)
        if identify_zeros && isz(a, model=model)
           coef = 0
        elseif symmetry == true
            if model == "Ising"
                a = reduce_mirror(a, N)
            else
                a = reduce_perm(a)
            end
        end
    end
    return a,coef
end

# binary search in a sorted sequence
function bfind(A, a)
    low = 1
    high = length(A)
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if A[mid] == a
           return mid
        elseif A[mid] < a
            low = mid + 1
        else
            high = mid - 1
        end
    end
    return nothing
end

function get_basis(N, d; label=1)
    if label == 1
        basis = [tuple(Int[], Vector{Int}[])]
        append!(basis, [tuple([3i-2], Vector{Int}[]) for i = 1:N])
        append!(basis, [tuple(Int[], [[3i-2]]) for i = 1:ceil(Int, N/2)])
        if d > 1
            for i = 1:N-1
                append!(basis, [tuple([3i-2;3j+1], Vector{Int}[]) for j = i:N-1])
                append!(basis, [tuple([3i-1;3j+2], Vector{Int}[]) for j = i:N-1])
                append!(basis, [tuple([3i;3j+3], Vector{Int}[]) for j = i:N-1])
                push!(basis, tuple([3i+1], [reduce_mirror([3i-2], N)]), tuple([3i-2], [reduce_mirror([3i+1], N)]))
            end
            # for i = 1:floor(Int, N/2)
            #     push!(basis, tuple(Int[], [[3i-2], [3i+1]]), tuple(Int[], [[3i-2;3i+1]]), tuple(Int[], [[3i-1;3i+2]]), tuple(Int[], [[3i;3i+3]]))
            # end
        end
        if d > 2
            for i = 1:N-2
               push!(basis, tuple([3i-2;3i+1;3i+4], Vector{Int}[]), tuple([3i-2;3i+2;3i+5], Vector{Int}[]), tuple([3i-1;3i+1;3i+5], Vector{Int}[]),
               tuple([3i-1;3i+2;3i+4], Vector{Int}[]), tuple([3i-2;3i+3;3i+6], Vector{Int}[]), tuple([3i;3i+1;3i+6], Vector{Int}[]), tuple([3i;3i+3;3i+4], Vector{Int}[]))
            end
        end
        if d > 3
            for i = 1:N-3
                push!(basis, tuple([3i-2;3i+1;3i+4;3i+7], Vector{Int}[]), tuple([3i-1;3i+2;3i+5;3i+8], Vector{Int}[]), tuple([3i;3i+3;3i+6;3i+9], Vector{Int}[]),
                tuple([3i-2;3i+1;3i+5;3i+8], Vector{Int}[]), tuple([3i-2;3i+2;3i+4;3i+8], Vector{Int}[]), tuple([3i-1;3i+1;3i+4;3i+8], Vector{Int}[]),
                tuple([3i-2;3i+2;3i+5;3i+7], Vector{Int}[]), tuple([3i-1;3i+1;3i+5;3i+7], Vector{Int}[]), tuple([3i-1;3i+2;3i+4;3i+7], Vector{Int}[]),
                tuple([3i-1;3i+2;3i+6;3i+9], Vector{Int}[]), tuple([3i-1;3i+3;3i+5;3i+9], Vector{Int}[]), tuple([3i;3i+2;3i+5;3i+9], Vector{Int}[]),
                tuple([3i-1;3i+3;3i+6;3i+8], Vector{Int}[]), tuple([3i;3i+2;3i+6;3i+8], Vector{Int}[]), tuple([3i;3i+3;3i+5;3i+8], Vector{Int}[]),
                tuple([3i;3i+3;3i+4;3i+7], Vector{Int}[]), tuple([3i;3i+1;3i+6;3i+7], Vector{Int}[]), tuple([3i-2;3i+3;3i+6;3i+7], Vector{Int}[]),
                tuple([3i;3i+1;3i+4;3i+9], Vector{Int}[]), tuple([3i-2;3i+3;3i+4;3i+9], Vector{Int}[]), tuple([3i-2;3i+1;3i+6;3i+9], Vector{Int}[]),)
            end
        end
        if d > 1
            for i = 1:N-1
                append!(basis, [tuple([3i-1;3j+3], Vector{Int}[]) for j = i:N-1])
                append!(basis, [tuple([3i;3j+2], Vector{Int}[]) for j = i:N-1])
            end
            # for i = 1:floor(Int, N/2)
            #     push!(basis, tuple(Int[], [[3i-1;3i+3]]), tuple(Int[], [[3i;3i+2]]))
            # end
        end
        if d > 2
            for i = 1:N-2
               push!(basis, tuple([3i-2;3i+2;3i+6], Vector{Int}[]), tuple([3i-1;3i+1;3i+6], Vector{Int}[]), tuple([3i-1;3i+3;3i+4], Vector{Int}[]), 
               tuple([3i-2;3i+3;3i+5], Vector{Int}[]), tuple([3i;3i+1;3i+5], Vector{Int}[]), tuple([3i;3i+2;3i+4], Vector{Int}[]))
            end
        end
        if d > 3
            for i = 1:N-3
                push!(basis, tuple([3i-2;3i+1;3i+5;3i+9], Vector{Int}[]), tuple([3i-2;3i+2;3i+4;3i+9], Vector{Int}[]), tuple([3i-1;3i+1;3i+4;3i+9], Vector{Int}[]),
                tuple([3i-2;3i+2;3i+6;3i+7], Vector{Int}[]), tuple([3i-1;3i+1;3i+6;3i+7], Vector{Int}[]), tuple([3i-1;3i+3;3i+4;3i+7], Vector{Int}[]),
                tuple([3i-2;3i+1;3i+6;3i+8], Vector{Int}[]), tuple([3i-2;3i+3;3i+4;3i+8], Vector{Int}[]), tuple([3i;3i+1;3i+4;3i+8], Vector{Int}[]),
                tuple([3i-2;3i+3;3i+5;3i+7], Vector{Int}[]), tuple([3i;3i+1;3i+5;3i+7], Vector{Int}[]), tuple([3i;3i+2;3i+4;3i+7], Vector{Int}[]),
                tuple([3i-1;3i+2;3i+5;3i+9], Vector{Int}[]), tuple([3i-1;3i+2;3i+6;3i+8], Vector{Int}[]), tuple([3i-1;3i+3;3i+5;3i+8], Vector{Int}[]), tuple([3i-2;3i+2;3i+5;3i+8], Vector{Int}[]),
                tuple([3i;3i+3;3i+6;3i+8], Vector{Int}[]), tuple([3i;3i+3;3i+5;3i+9], Vector{Int}[]), tuple([3i;3i+2;3i+6;3i+9], Vector{Int}[]), tuple([3i-1;3i+3;3i+6;3i+9], Vector{Int}[]))
            end
        end
    else
        basis = [tuple([3i-1], Vector{Int}[]) for i = 1:N]
        if d > 1
            append!(basis, [tuple([3i-2;3i+2], Vector{Int}[])  for i = 1:N-1])
            append!(basis, [tuple([3i-1;3i+1], Vector{Int}[]) for i = 1:N-1])
            # for i = 1:N-1
            #     push!(basis, tuple([3i+2], [reduce_mirror([3i-2], N)]), tuple([3i-1], [reduce_mirror([3i+1], N)]))
            # end
        end
        if d > 2
            for i = 1:N-2
               push!(basis, tuple([3i-1;3i+2;3i+5], Vector{Int}[]), tuple([3i-1;3i+1;3i+4], Vector{Int}[]), tuple([3i-2;3i+2;3i+4], Vector{Int}[]), tuple([3i-2;3i+1;3i+5], Vector{Int}[]),
               tuple([3i-1;3i+3;3i+6], Vector{Int}[]), tuple([3i;3i+2;3i+6], Vector{Int}[]), tuple([3i;3i+3;3i+5], Vector{Int}[]))
            end
        end
        if d > 3
            for i = 1:N-3
                push!(basis, tuple([3i;3i+3;3i+4;3i+8], Vector{Int}[]), tuple([3i;3i+1;3i+6;3i+8], Vector{Int}[]), tuple([3i-2;3i+3;3i+6;3i+8], Vector{Int}[]),
                tuple([3i-2;3i+3;3i+5;3i+9], Vector{Int}[]), tuple([3i;3i+1;3i+5;3i+9], Vector{Int}[]), tuple([3i-2;3i+2;3i+6;3i+9], Vector{Int}[]),
                tuple([3i;3i+3;3i+5;3i+7], Vector{Int}[]), tuple([3i;3i+2;3i+6;3i+7], Vector{Int}[]), tuple([3i-1;3i+3;3i+6;3i+7], Vector{Int}[]),
                tuple([3i-1;3i+3;3i+4;3i+9], Vector{Int}[]), tuple([3i;3i+2;3i+4;3i+9], Vector{Int}[]), tuple([3i-1;3i+1;3i+6;3i+9], Vector{Int}[]),
                tuple([3i-2;3i+1;3i+4;3i+8], Vector{Int}[]), tuple([3i-2;3i+1;3i+5;3i+7], Vector{Int}[]), tuple([3i-2;3i+2;3i+4;3i+7], Vector{Int}[]), tuple([3i-1;3i+1;3i+4;3i+7], Vector{Int}[]),
                tuple([3i-1;3i+2;3i+5;3i+7], Vector{Int}[]), tuple([3i-1;3i+2;3i+4;3i+8], Vector{Int}[]), tuple([3i-1;3i+1;3i+5;3i+8], Vector{Int}[]), tuple([3i-2;3i+2;3i+5;3i+8], Vector{Int}[]))
            end
        end
        append!(basis, [tuple([3i], Vector{Int}[]) for i = 1:N])
        if d > 1
            append!(basis, [tuple([3i-2;3i+3], Vector{Int}[]) for i = 1:N-1])
            append!(basis, [tuple([3i;3i+1], Vector{Int}[]) for i = 1:N-1])
            # for i = 1:N-1
            #     push!(basis, tuple([3i+3], [reduce_mirror([3i-2], N)]), tuple([3i], [reduce_mirror([3i+1], N)]))
            # end
        end
        if d > 2
            for i = 1:N-2
                push!(basis, tuple([3i;3i+3;3i+6], Vector{Int}[]), tuple([3i;3i+1;3i+4], Vector{Int}[]), tuple([3i-2;3i+3;3i+4], Vector{Int}[]), 
                tuple([3i-2;3i+1;3i+6], Vector{Int}[]), tuple([3i;3i+2;3i+5], Vector{Int}[]), tuple([3i-1;3i+3;3i+5], Vector{Int}[]), tuple([3i-1;3i+2;3i+6], Vector{Int}[]))
            end
        end
        if d > 3
            for i = 1:N-3
                push!(basis, tuple([3i-1;3i+2;3i+4;3i+9], Vector{Int}[]), tuple([3i-1;3i+1;3i+5;3i+9], Vector{Int}[]), tuple([3i-2;3i+2;3i+5;3i+9], Vector{Int}[]),
                tuple([3i-1;3i+1;3i+6;3i+8], Vector{Int}[]), tuple([3i-2;3i+2;3i+6;3i+8], Vector{Int}[]), tuple([3i-2;3i+3;3i+5;3i+8], Vector{Int}[]),
                tuple([3i-1;3i+2;3i+6;3i+7], Vector{Int}[]), tuple([3i-1;3i+3;3i+5;3i+7], Vector{Int}[]), tuple([3i;3i+2;3i+5;3i+7], Vector{Int}[]),
                tuple([3i;3i+2;3i+4;3i+8], Vector{Int}[]), tuple([3i-1;3i+3;3i+4;3i+8], Vector{Int}[]), tuple([3i;3i+1;3i+5;3i+8], Vector{Int}[]),
                tuple([3i-2;3i+1;3i+4;3i+9], Vector{Int}[]), tuple([3i-2;3i+1;3i+6;3i+7], Vector{Int}[]), tuple([3i-2;3i+3;3i+4;3i+7], Vector{Int}[]), tuple([3i;3i+1;3i+4;3i+7], Vector{Int}[]),             
                tuple([3i;3i+3;3i+6;3i+7], Vector{Int}[]), tuple([3i;3i+3;3i+4;3i+9], Vector{Int}[]), tuple([3i;3i+1;3i+6;3i+9], Vector{Int}[]), tuple([3i-2;3i+3;3i+6;3i+9], Vector{Int}[]))
            end
        end
    end
    return basis
end

function get_bulkbasis(N, d; label=1)
    if label == 1
        basis = [tuple([3i-2], Vector{Int}[]) for i = 2:N-1]
        append!(basis, [tuple(Int[], [[3i-2]]) for i = 2:ceil(Int, N/2)])
        if d > 1
            for i = 2:N-2
                append!(basis, [tuple([3i-2;3j+1], Vector{Int}[]) for j = i:N-2])
                append!(basis, [tuple([3i-1;3j+2], Vector{Int}[]) for j = i:N-2])
                append!(basis, [tuple([3i;3j+3], Vector{Int}[]) for j = i:N-2])
                push!(basis, tuple([3i+1], [reduce_mirror([3i-2], N)]), tuple([3i-2], [reduce_mirror([3i+1], N)]))
            end
            # for i = 2:floor(Int, N/2)
            #     push!(basis, tuple(Int[], [[3i-2], [3i+1]]), tuple(Int[], [[3i-2;3i+1]]), tuple(Int[], [[3i-1;3i+2]]), tuple(Int[], [[3i;3i+3]]))
            # end
        end
        if d > 2
            for i = 2:N-3
               push!(basis, tuple([3i-2;3i+1;3i+4], Vector{Int}[]), tuple([3i-2;3i+2;3i+5], Vector{Int}[]), tuple([3i-1;3i+1;3i+5], Vector{Int}[]),
               tuple([3i-1;3i+2;3i+4], Vector{Int}[]), tuple([3i-2;3i+3;3i+6], Vector{Int}[]), tuple([3i;3i+1;3i+6], Vector{Int}[]), tuple([3i;3i+3;3i+4], Vector{Int}[]))
            end
        end
        if d > 1
            for i = 2:N-2
                append!(basis, [tuple([3i-1;3j+3], Vector{Int}[]) for j = i:N-2])
                append!(basis, [tuple([3i;3j+2], Vector{Int}[]) for j = i:N-2])
            end
            # for i = 2:floor(Int, N/2)
            #     push!(basis, tuple(Int[], [[3i-1;3i+3]]), tuple(Int[], [[3i;3i+2]]))
            # end
        end
        if d > 2
            for i = 2:N-3
               push!(basis, tuple([3i-2;3i+2;3i+6], Vector{Int}[]), tuple([3i-1;3i+1;3i+6], Vector{Int}[]), tuple([3i-1;3i+3;3i+4], Vector{Int}[]), 
               tuple([3i-2;3i+3;3i+5], Vector{Int}[]), tuple([3i;3i+1;3i+5], Vector{Int}[]), tuple([3i;3i+2;3i+4], Vector{Int}[]))
            end
        end
    else
        basis = [tuple([3i-1], Vector{Int}[]) for i = 2:N-1]
        if d > 1
            append!(basis, [tuple([3i-2;3i+2], Vector{Int}[])  for i = 2:N-2])
            append!(basis, [tuple([3i-1;3i+1], Vector{Int}[]) for i = 2:N-2])
            # for i = 2:N-2
            #     push!(basis, tuple([3i+2], [reduce_mirror([3i-2], N)]), tuple([3i-1], [reduce_mirror([3i+1], N)]))
            # end
        end
        if d > 2
            for i = 2:N-3
               push!(basis, tuple([3i-1;3i+2;3i+5], Vector{Int}[]), tuple([3i-1;3i+1;3i+4], Vector{Int}[]), tuple([3i-2;3i+2;3i+4], Vector{Int}[]), tuple([3i-2;3i+1;3i+5], Vector{Int}[]),
               tuple([3i-1;3i+3;3i+6], Vector{Int}[]), tuple([3i;3i+2;3i+6], Vector{Int}[]), tuple([3i;3i+3;3i+5], Vector{Int}[]))
            end
        end
        append!(basis, [tuple([3i], Vector{Int}[]) for i = 2:N-1])
        if d > 1
            append!(basis, [tuple([3i-2;3i+3], Vector{Int}[]) for i = 2:N-2])
            append!(basis, [tuple([3i;3i+1], Vector{Int}[]) for i = 2:N-2])
            # for i = 2:N-2
            #     push!(basis, tuple([3i+3], [reduce_mirror([3i-2], N)]), tuple([3i], [reduce_mirror([3i+1], N)]))
            # end
        end
        if d > 2
            for i = 2:N-3
                push!(basis, tuple([3i;3i+3;3i+6], Vector{Int}[]), tuple([3i;3i+1;3i+4], Vector{Int}[]), tuple([3i-2;3i+3;3i+4], Vector{Int}[]), 
                tuple([3i-2;3i+1;3i+6], Vector{Int}[]), tuple([3i;3i+2;3i+5], Vector{Int}[]), tuple([3i-1;3i+3;3i+5], Vector{Int}[]), tuple([3i-1;3i+2;3i+6], Vector{Int}[]))
            end
        end
    end
    return basis
end

function get_kagome_basis(N, triples, edges, d; label=1)
    if label == 1
        basis = [tuple(Int[], Vector{Int}[])]
        for i = 1:N-1
            append!(basis, [tuple([3i-2;3j+1], Vector{Int}[]) for j = i:N-1])
            append!(basis, [tuple([3i-1;3j+2], Vector{Int}[]) for j = i:N-1])
            append!(basis, [tuple([3i;3j+3], Vector{Int}[]) for j = i:N-1])
        end
        if d > 2
            for a in triples
                push!(basis, tuple([3*a[1]-2;3*a[2]-1;3*a[3]], Vector{Int}[]), tuple([3*a[1]-2;3*a[2];3*a[3]-1], Vector{Int}[]), tuple([3*a[1]-1;3*a[2]-2;3*a[3]], Vector{Int}[]),
                tuple([3*a[1]-1;3*a[2];3*a[3]-2], Vector{Int}[]), tuple([3*a[1];3*a[2]-2;3*a[3]-1], Vector{Int}[]), tuple([3*a[1];3*a[2]-1;3*a[3]-2], Vector{Int}[]))
            end
        end
    elseif label == 2
        basis = [tuple([3i-2], Vector{Int}[]) for i = 1:N]
        append!(basis, [tuple(Int[], [[3i-2]]) for i = 1:N])
        if d > 2
            for a in triples
                push!(basis, tuple([3*a[1]-2;3*a[2]-2;3*a[3]-2], Vector{Int}[]), tuple([3*a[1]-2;3*a[2]-1;3*a[3]-1], Vector{Int}[]), tuple([3*a[1]-1;3*a[2]-2;3*a[3]-1], Vector{Int}[]),
                tuple([3*a[1]-1;3*a[2]-1;3*a[3]-2], Vector{Int}[]), tuple([3*a[1]-2;3*a[2];3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2]-2;3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2];3*a[3]-2], Vector{Int}[]))
            end
        end
        for a in edges
            append!(basis, [tuple([3*a[1]-1;3*a[2]], Vector{Int}[]), tuple([3*a[1];3*a[2]-1], Vector{Int}[])])
        end
        for a in triples
            append!(basis, [tuple([3*a[1]-1;3*a[2]], Vector{Int}[]), tuple([3*a[1]-1;3*a[3]], Vector{Int}[]), tuple([3*a[2]-1;3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2]-1], Vector{Int}[]), tuple([3*a[1];3*a[3]-1], Vector{Int}[]), tuple([3*a[2];3*a[3]-1], Vector{Int}[])])
        end
    elseif label == 3
        basis = [tuple([3i-1], Vector{Int}[]) for i = 1:N]
        append!(basis, [tuple(Int[], [[3i-1]]) for i = 1:N])
        if d > 2
            for a in triples
                push!(basis, tuple([3*a[1]-1;3*a[2]-1;3*a[3]-1], Vector{Int}[]), tuple([3*a[1]-1;3*a[2]-2;3*a[3]-2], Vector{Int}[]), tuple([3*a[1]-2;3*a[2]-1;3*a[3]-2], Vector{Int}[]),
                tuple([3*a[1]-2;3*a[2]-2;3*a[3]-1], Vector{Int}[]), tuple([3*a[1]-1;3*a[2];3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2]-1;3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2];3*a[3]-1], Vector{Int}[]))
            end
        end
        for a in edges
            append!(basis, [tuple([3*a[1]-2;3*a[2]], Vector{Int}[]), tuple([3*a[1];3*a[2]-2], Vector{Int}[])])
        end
        for a in triples
            append!(basis, [tuple([3*a[1]-2;3*a[2]], Vector{Int}[]), tuple([3*a[1]-2;3*a[3]], Vector{Int}[]), tuple([3*a[2]-2;3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2]-2], Vector{Int}[]), tuple([3*a[1];3*a[3]-2], Vector{Int}[]), tuple([3*a[2];3*a[3]-2], Vector{Int}[])])
        end
    else
        basis = [tuple([3i], Vector{Int}[]) for i = 1:N]
        append!(basis, [tuple(Int[], [[3i]]) for i = 1:N])
        if d > 2
            for a in triples
                push!(basis, tuple([3*a[1];3*a[2];3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2]-2;3*a[3]-2], Vector{Int}[]), tuple([3*a[1]-2;3*a[2];3*a[3]-2], Vector{Int}[]),
                tuple([3*a[1]-2;3*a[2]-2;3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2]-1;3*a[3]-1], Vector{Int}[]), tuple([3*a[1]-1;3*a[2];3*a[3]-1], Vector{Int}[]), tuple([3*a[1]-1;3*a[2]-1;3*a[3]], Vector{Int}[]))
            end
        end
        for a in edges
            append!(basis, [tuple([3*a[1]-1;3*a[2]-2], Vector{Int}[]), tuple([3*a[1]-2;3*a[2]-1], Vector{Int}[])])
        end
        for a in triples
            append!(basis, [tuple([3*a[1]-1;3*a[2]-2], Vector{Int}[]), tuple([3*a[1]-1;3*a[3]-2], Vector{Int}[]), tuple([3*a[2]-1;3*a[3]-2], Vector{Int}[]), tuple([3*a[1]-2;3*a[2]-1], Vector{Int}[]), tuple([3*a[1]-2;3*a[3]-1], Vector{Int}[]), tuple([3*a[2]-2;3*a[3]-1], Vector{Int}[])])
        end
    end
    return basis
end

function get_kagome_bulkbasis(N, triples, edges, d; label=1)
    if label == 1
        basis = Tuple{Vector{Int}, Vector{Vector{Int}}}[]
        if d > 1
            for a in edges
                append!(basis, [tuple([3*a[1]-2;3*a[2]-2], Vector{Int}[]), tuple([3*a[1]-1;3*a[2]-1], Vector{Int}[]), tuple([3*a[1];3*a[2]], Vector{Int}[])])
            end
            for a in triples
                append!(basis, [tuple([3*a[1]-2;3*a[2]-2], Vector{Int}[]), tuple([3*a[1]-1;3*a[2]-1], Vector{Int}[]), tuple([3*a[1];3*a[2]], Vector{Int}[]), tuple([3*a[1]-2;3*a[3]-2], Vector{Int}[]), 
                tuple([3*a[1]-1;3*a[3]-1], Vector{Int}[]), tuple([3*a[1];3*a[3]], Vector{Int}[]), tuple([3*a[2]-2;3*a[3]-2], Vector{Int}[]), tuple([3*a[2]-1;3*a[3]-1], Vector{Int}[]), tuple([3*a[2];3*a[3]], Vector{Int}[])])
            end
        end
    elseif label == 2
        if N == 5
            basis = [tuple([1], Vector{Int}[])]
        elseif N == 13
            basis = [tuple([3i-2], Vector{Int}[]) for i = 1:5]
        elseif N == 27
            basis = [tuple([3i-2], Vector{Int}[]) for i = 1:13]
        elseif N == 45
            basis = [tuple([3i-2], Vector{Int}[]) for i = 1:27]
        else
            @error "Wrong number of sites!"
        end
        if d > 1
            for a in edges
                append!(basis, [tuple([3*a[1]-1;3*a[2]], Vector{Int}[]), tuple([3*a[1];3*a[2]-1], Vector{Int}[])])
            end
            for a in triples
                append!(basis, [tuple([3*a[1]-1;3*a[2]], Vector{Int}[]), tuple([3*a[1]-1;3*a[3]], Vector{Int}[]), tuple([3*a[2]-1;3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2]-1], Vector{Int}[]), tuple([3*a[1];3*a[3]-1], Vector{Int}[]), tuple([3*a[2];3*a[3]-1], Vector{Int}[])])
            end
        end
    elseif label == 3
        if N == 5
            basis = [tuple([2], Vector{Int}[])]
        elseif N == 13
            basis = [tuple([3i-1], Vector{Int}[]) for i = 1:5]
        elseif N == 27
            basis = [tuple([3i-1], Vector{Int}[]) for i = 1:13]
        elseif N == 45
            basis = [tuple([3i-1], Vector{Int}[]) for i = 1:27]
        else
            @error "Wrong number of sites!"
        end
        if d > 1
            for a in edges
                append!(basis, [tuple([3*a[1]-2;3*a[2]], Vector{Int}[]), tuple([3*a[1];3*a[2]-2], Vector{Int}[])])
            end
            for a in triples
                append!(basis, [tuple([3*a[1]-2;3*a[2]], Vector{Int}[]), tuple([3*a[1]-2;3*a[3]], Vector{Int}[]), tuple([3*a[2]-2;3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2]-2], Vector{Int}[]), tuple([3*a[1];3*a[3]-2], Vector{Int}[]), tuple([3*a[2];3*a[3]-2], Vector{Int}[])])
            end
        end
    else
        if N == 5
            basis = [tuple([3], Vector{Int}[])]
        elseif N == 13
            basis = [tuple([3i], Vector{Int}[]) for i = 1:5]
        elseif N == 27
            basis = [tuple([3i], Vector{Int}[]) for i = 1:13]
        elseif N == 45
            basis = [tuple([3i], Vector{Int}[]) for i = 1:27]
        else
            @error "Wrong number of sites!"
        end
        if d > 1
            for a in edges
                append!(basis, [tuple([3*a[1]-2;3*a[2]-1], Vector{Int}[]), tuple([3*a[1]-1;3*a[2]-2], Vector{Int}[])])
            end
            for a in triples
                append!(basis, [tuple([3*a[1]-2;3*a[2]-1], Vector{Int}[]), tuple([3*a[1]-2;3*a[3]-1], Vector{Int}[]), tuple([3*a[2]-2;3*a[3]-1], Vector{Int}[]), tuple([3*a[1]-1;3*a[2]-2], Vector{Int}[]), tuple([3*a[1]-1;3*a[3]-2], Vector{Int}[]), tuple([3*a[2]-1;3*a[3]-2], Vector{Int}[])])
            end
        end
    end
    return basis
end

function arrange(Locb, coef)
    nLocb = copy(Locb)
    sort!(nLocb)
    unique!(nLocb)
    ncoef = zeros(typeof(coef[1]), length(nLocb))
    for (i, item) in enumerate(Locb)
        ind = bfind(nLocb, item)
        ncoef[ind] += coef[i]
    end
    ind = ncoef .!= 0
    return nLocb[ind],ncoef[ind]
end

function update!(coef, bis, bi, N; λ=1, model="Ising", realify=false, identify_zeros=true, symmetry=true)
    bi,c = reduce!(bi, N, model=model, realify=realify, identify_zeros=identify_zeros, symmetry=symmetry)
    if c != 0
        push!(coef, λ*c)
        push!(bis, bi)
    end
end

function count_reduce(item, ba)
    loc = ceil.(Int, ba/3)
    u = bfind(loc, ceil(Int, item[1]/3))
    if length(item) == 1
        return u !== nothing && mod(ba[u], 3) != mod(item[1], 3)
    else
        v = bfind(loc, ceil(Int, item[2]/3))
        return (u !== nothing && mod(ba[u], 3) != mod(item[1], 3)) + (v !== nothing && mod(ba[v], 3) != mod(item[2], 3)) == 1
    end
end

function PSDstate_entry(ba1, ba2, H::ncpoly, L; model="Ising", realify=false, identify_zeros=true, symmetry=true)
    if realify == true
        coef = Float16[]
    else
        coef = ComplexF16[]
    end
    bis = Vector{Int}[]
    for (i, item) in enumerate(H.supp)
        s = count_reduce(item, ba1)
        t = count_reduce(item, ba2)
        if s + t > 0
            bi = [ba1; item; ba2]
            update!(coef, bis, bi, L, λ=H.coe[i]*(s+t), model=model, realify=realify, identify_zeros=identify_zeros, symmetry=symmetry)
        end
    end
    if !isempty(coef)
        bis, coef = arrange(bis, coef)
    end
    return bis, coef
end

function generate_mons(N, d)
    mons = Vector{Int}[]
    for i = 2:N-2, j = i+1:N-1, u = 1:3, v = 1:3
        push!(mons, [3*(i-1)+u; 3*(j-1)+v])
    end
    if d > 2
        for i = 2:N-3, j = i+1:N-2, k = j+1:N-1, u = 1:3, v = 1:3, s = 1:3
            push!(mons, [3*(i-1)+u; 3*(j-1)+v; 3*(k-1)+s])
        end
    end
    if d > 3
        for i = 2:N-4, j = i+1:N-3, k = j+1:N-2, l = k+1:N-1, u = 1:3, v = 1:3, s = 1:3, t = 1:3
            push!(mons, [3*(i-1)+u; 3*(j-1)+v; 3*(k-1)+s; 3*(l-1)+t])
        end
    end
    if d > 4
        for i = 2:N-5, j = i+1:N-4, k = j+1:N-3, l = k+1:N-2, p = l+1:N-1, u = 1:3, v = 1:3, s = 1:3, t = 1:3, q = 1:3
            push!(mons, [3*(i-1)+u; 3*(j-1)+v; 3*(k-1)+s; 3*(l-1)+t; 3*(p-1)+q])
        end
    end
    if d > 5
        for i = 2:N-6, j = i+1:N-5, k = j+1:N-4, l = k+1:N-3, p = l+1:N-2, a = p+1:N-1, u = 1:3, v = 1:3, s = 1:3, t = 1:3, q = 1:3, b = 1:3
            push!(mons, [3*(i-1)+u; 3*(j-1)+v; 3*(k-1)+s; 3*(l-1)+t; 3*(p-1)+q; 3*(a-1)+b])
        end
    end
    if d > 6
        for i = 2:N-7, j = i+1:N-6, k = j+1:N-5, l = k+1:N-4, p = l+1:N-3, a = p+1:N-2, c = a+1:N-1, u = 1:3, v = 1:3, s = 1:3, t = 1:3, q = 1:3, b = 1:3, d = 1:3
            push!(mons, [3*(i-1)+u; 3*(j-1)+v; 3*(k-1)+s; 3*(l-1)+t; 3*(p-1)+q; 3*(a-1)+b; 3*(c-1)+d])
        end
    end
    return mons
end

function generate_mons_kagome(N, d)
    mons = Vector{Int}[]
    for i = 1:N-1, j = i+1:N, u = 1:3, v = 1:3
        push!(mons, [3*(i-1)+u; 3*(j-1)+v])
    end
    if d > 2
        for i = 1:N-2, j = i+1:N-1, k = j+1:N, u = 1:3, v = 1:3, s = 1:3
            push!(mons, [3*(i-1)+u; 3*(j-1)+v; 3*(k-1)+s])
        end
    end
    if d > 3
        for i = 1:N-3, j = i+1:N-2, k = j+1:N-1, l = k+1:N, u = 1:3, v = 1:3, s = 1:3, t = 1:3
            push!(mons, [3*(i-1)+u; 3*(j-1)+v; 3*(k-1)+s; 3*(l-1)+t])
        end
    end
    if d > 4
        for i = 1:N-4, j = i+1:N-3, k = j+1:N-2, l = k+1:N-1, p = l+1:N, u = 1:3, v = 1:3, s = 1:3, t = 1:3, q = 1:3
            push!(mons, [3*(i-1)+u; 3*(j-1)+v; 3*(k-1)+s; 3*(l-1)+t; 3*(p-1)+q])
        end
    end
    if d > 5
        for i = 1:N-5, j = i+1:N-4, k = j+1:N-3, l = k+1:N-2, p = l+1:N-1, a = p+1:N, u = 1:3, v = 1:3, s = 1:3, t = 1:3, q = 1:3, b = 1:3
            push!(mons, [3*(i-1)+u; 3*(j-1)+v; 3*(k-1)+s; 3*(l-1)+t; 3*(p-1)+q; 3*(a-1)+b])
        end
    end
    if d > 6
        for i = 1:N-6, j = i+1:N-5, k = j+1:N-4, l = k+1:N-3, p = l+1:N-2, a = p+1:N-1, c = a+1:N, u = 1:3, v = 1:3, s = 1:3, t = 1:3, q = 1:3, b = 1:3, d = 1:3
            push!(mons, [3*(i-1)+u; 3*(j-1)+v; 3*(k-1)+s; 3*(l-1)+t; 3*(p-1)+q; 3*(a-1)+b; 3*(c-1)+d])
        end
    end
    return mons
end

function filter_mons(mons, tsupp, H::ncpoly, L; model="Ising", identify_zeros=true, symmetry=true)
    lc = zeros(Float64, length(mons))
    rd = ceil.(Int, rand(length(tsupp))*10^8)
    for (k, mon) in enumerate(mons)
        flag = 1
        for (i, item) in enumerate(H.supp)
            bi,c = reduce!([item; mon], L, model=model, identify_zeros=identify_zeros, symmetry=symmetry)
            if imag(c) != 0
                Locb = bfind(tsupp, [bi])
                if Locb === nothing
                    flag = 0
                    break
                end
                lc[k] += H.coe[i]*imag(c)*rd[Locb]
            end
        end
        if flag == 0
            lc[k] = 0
        end
    end
    ind = unique(i -> abs(lc[i]), eachindex(lc))
    return mons[ind[abs.(lc[ind]) .> 0]]
end
