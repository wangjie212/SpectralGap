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
function isz(a::Vector{Int})
    return any(i->isodd(count(isequal(i), mod.(a,3))), [0, 2])
end

# reduction w.r.t mirror symmetry
function reduce4(a::Vector{Int}, N)
    ra = reverse(a)
    ma = 3*(N .- ceil.(Int, ra/3)) + smod.(ra, 3)
    return min(a, ma)
end

function smod(i, s)
    r = mod(i, s)
    return r == 0 ? s : r
end

# implement all reductions
function reduce!(a::Vector{Int}, N; realify=false)
    reduce1!(a)
    reduce3!(a)
    a,coef = reduce2!(a, realify=realify)
    reduce3!(a)
    if !isempty(a)
        if isz(a)
           coef = 0
        else
           a = reduce4(a, N)
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

function get_wbasis(N, d; label=1)
    if label == 1
        basis = [[3*i-2] for i = 2:N-1]
        if d > 1
            append!(basis, [[3*i-2;3*i+1] for i = 2:N-2])
            append!(basis, [[3*i-1;3*i+2] for i = 2:N-2])
            append!(basis, [[3*i;3*i+3] for i = 2:N-2])
        end
    else
        # basis = [[[3*i-1] for i = 2:N-1]; [[3*i] for i = 2:N-1]]
        basis = [[3*i-1] for i = 2:N-1]
        if d > 1
            append!(basis, [[3*i-2;3*i+2] for i = 2:N-2])
            append!(basis, [[3*i-1;3*i+1] for i = 2:N-2])
        end
        append!(basis, [[3*i] for i = 2:N-1])
        if d > 1
            append!(basis, [[3*i-2;3*i+3] for i = 2:N-2])
            append!(basis, [[3*i;3*i+1] for i = 2:N-2])
        end
    end
    return basis
end

function get_sbasis(N)
    basis = [[Vector{Int}[]]; [[[3*i-2]] for i = 1:ceil(Int, N/2)]]
    return basis
end

function get_basis(N, d; label=1)
    if label == 1
        basis = [tuple(Int[], Vector{Int}[])]
        append!(basis, [tuple([3*i-2], Vector{Int}[]) for i = 1:N])
        append!(basis, [tuple(Int[], [[3*i-2]]) for i = 1:ceil(Int, N/2)])
        append!(basis, [tuple([3*i-2;3*i+1], Vector{Int}[]) for i = 1:N-1])
        for i = 1:N-1
            append!(basis, [tuple([3*i-1;3*j+2], Vector{Int}[]) for j = i:N-1])
            append!(basis, [tuple([3*i;3*j+3], Vector{Int}[]) for j = i:N-1])
            # append!(basis, [tuple([3*i-2], [[3*j+1]]) for j = i:N-1])
            # append!(basis, [tuple([3*i+1], [[3*j-2]]) for j = i:N-1])
            # append!(basis, [tuple(Int[], [[3*i-2;3*j+1]]) for j = i:N-1])
            # append!(basis, [tuple(Int[], [[3*i-1;3*j+2]]) for j = i:N-1])
            # append!(basis, [tuple(Int[], [[3*i;3*j+3]]) for j = i:N-1])
            # append!(basis, [tuple(Int[], [[3*i-2], [3*j+1]]) for j = i:N-1])
            # append!(basis, [tuple([3*i-1;3*j+3], Vector{Int}[]) for j = i:N-1])
            # append!(basis, [tuple([3*i;3*j+2], Vector{Int}[]) for j = i:N-1])
            # push!(basis, tuple(Int[], [[3*i-1;3*i+2]]), tuple(Int[], [[3*i;3*i+3]]), tuple([3*i+1], [[3*i-2]]), 
            # tuple([3*i-2], [[3*i+1]]), tuple(Int[], [[3*i-2;3*i+1]]), tuple(Int[], [[3*i-2], [3*i+1]]))
        end
        # for i = 1:N-1
        #     push!(basis, tuple([3*i-1;3*i+3], Vector{Int}[]), tuple([3*i;3*i+2], Vector{Int}[]))
        # end
        if d > 2
            append!(basis, [tuple(Int[], [[3*i-2;3*i+1]]) for i = 1:floor(Int, N/2)])
            append!(basis, [tuple(Int[], [[3*i;3*i+3]]) for i = 1:floor(Int, N/2)])
            for i = 1:floor(Int, N/2)
                append!(basis, [tuple(Int[], [[3*i-1;3*j+2]]) for j = i:N-1])
            end
            for i = 1:ceil(Int, N/2)
                append!(basis, [tuple(Int[], [[3*i-2], [3*j-2]]) for j = i:ceil(Int, N/2)])
            end
            for i = 1:N-1
                append!(basis, [tuple([3*i-2;3*j+1], Vector{Int}[]) for j = i:N-1])
            end
            for i = 1:N-2
               push!(basis, tuple([3*i-2;3*i+1;3*i+4], Vector{Int}[]), tuple([3*i-2;3*i+2;3*i+5], Vector{Int}[]), tuple([3*i-1;3*i+1;3*i+5], Vector{Int}[]),
               tuple([3*i-1;3*i+2;3*i+4], Vector{Int}[]), tuple([3*i-2;3*i+3;3*i+6], Vector{Int}[]), tuple([3*i;3*i+1;3*i+6], Vector{Int}[]), tuple([3*i;3*i+3;3*i+4], Vector{Int}[]))
            end
    #     end
    # elseif label == 2
    #     basis = [tuple(Int[], [[3*i-2]]) for i = 1:ceil(Int, N/2)]
    #     if d > 2
            for i = 1:floor(Int, N/2)
                append!(basis, [tuple(Int[], [[3*i;3*j]]) for j = i+1:N])
            end
            for i = 1:N-2, j = i+1:N-1, k = j+1:N
                push!(basis, tuple([3*i-2;3*j;3*k], Vector{Int}[]), tuple([3*i;3*j-2;3*k], Vector{Int}[]), tuple([3*i;3*j;3*k-2], Vector{Int}[]))
            end
            for i = 1:N-1
                push!(basis, tuple(Int[], [[3*i-2;3*i+1]]), tuple(Int[], [[3*i-1;3*i+2]]), tuple(Int[], [[3*i;3*i+3]]))
            end
            # for i = 1:N-2
            #    push!(basis, tuple(Int[3*i+1;3*i+4], [[3*i-2]]), tuple(Int[3*i+2;3*i+5], [[3*i-2]]), tuple(Int[3*i+3;3*i+6], [[3*i-2]]),
            #    tuple(Int[3*i-2;3*i+1], [[3*i+4]]), tuple(Int[3*i-1;3*i+2], [[3*i+4]]), tuple(Int[3*i;3*i+3], [[3*i+4]]))
            # end
        end
    else
        basis = [tuple([3*i-1], Vector{Int}[]) for i = 1:N]
        for i = 1:N-1
            push!(basis, tuple([3*i-2;3*i+2], Vector{Int}[]), tuple([3*i-1], [[3*i+1]]), tuple([3*i-2], [[3*i+2]]), tuple([3*i-1;3*i+1], Vector{Int}[]))
        end
        # basis = [tuple([3*i], Vector{Int}[]) for i = 1:N]
        append!(basis, [tuple([3*i], Vector{Int}[]) for i = 1:N])
        append!(basis, [tuple([3*i-2;3*i+3], Vector{Int}[]) for i = 1:N-1])
        append!(basis, [tuple([3*i;3*i+1], Vector{Int}[]) for i = 1:N-1])
        if d > 2
            for i = 1:N-2
               push!(basis, tuple([3*i;3*i+3;3*i+6], Vector{Int}[]), tuple([3*i;3*i+1;3*i+4], Vector{Int}[]), tuple([3*i-2;3*i+3;3*i+4], Vector{Int}[]),
               tuple([3*i-2;3*i+1;3*i+6], Vector{Int}[]), tuple([3*i;3*i+2;3*i+5], Vector{Int}[]), tuple([3*i-1;3*i+3;3*i+5], Vector{Int}[]), tuple([3*i-1;3*i+2;3*i+6], Vector{Int}[]))
            end
        end
    end
    return basis
end
