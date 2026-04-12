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
# function isz(a::Vector{Int})
#     return any(i->isodd(count(isequal(i), mod.(a,3))), [0, 2])
# end

function isz(a::Vector{Int})
    return any(i->isodd(count(isequal(i), mod.(a,3))), [0, 1, 2])
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
        # else
        #    a = reduce4(a, N)
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
        basis = [[3i-2] for i = 2:N-1]
    else
        basis = [[[3i-1] for i = 2:N-1]; [[3i] for i = 2:N-1]]
        if d > 1
            append!(basis, [[3i-2;3i+3] for i = 2:N-2])
            append!(basis, [[3i;3i+1] for i = 2:N-2])
        end
    end
    return basis
end

function get_sbasis(N)
    basis = [[Vector{Int}[]]; [[[3i-2]] for i = 1:ceil(Int, N/2)]]
    return basis
end

function get_basis(N, d; label=1)
    if label == 1
        basis = [tuple(Int[], Vector{Int}[])]
        append!(basis, [tuple([3i-2], Vector{Int}[]) for i = 2:N-1])
        for i = 2:N-2
            append!(basis, [tuple([3i-1;3j+2], Vector{Int}[]) for j = i:N-2])
        end
        for i = 1:N-1
            append!(basis, [tuple([3i;3j+3], Vector{Int}[]) for j = i:N-1])
        end
        if d > 2
            for i = 2:N-2
                append!(basis, [tuple([3i-2;3j+1], Vector{Int}[]) for j = i:N-2])
            end
        end
        append!(basis, [tuple([3i;3i+2], Vector{Int}[]) for i = 1:N-2])
        append!(basis, [tuple([3i-1;3i+3], Vector{Int}[]) for i = 2:N-1])
    else
        basis = [tuple(Int[1], Vector{Int}[]), tuple(Int[3N-2], Vector{Int}[]), tuple(Int[], [[1]])]
        append!(basis, [tuple([2;3i-1], Vector{Int}[]) for i = 2:N-1])
        append!(basis, [tuple([3i-1;3N-1], Vector{Int}[]) for i = 2:N-1])
        if d > 2
            append!(basis, [tuple(Int[], [[1], [3i-2]]) for i = 2:ceil(Int, N/2)])
            append!(basis, [tuple([1;3i-2], Vector{Int}[]) for i = 2:N-1])
            append!(basis, [tuple([3i-2;3N-2], Vector{Int}[]) for i = 2:N-1])
        end
        push!(basis, tuple([2;6], Vector{Int}[]), tuple([3N-3;3N-1], Vector{Int}[]))
    end
    return basis
end

function get_kagome_basis(N, triples, edges, d; label=1)
    if label == 1
        basis = [tuple(Int[], Vector{Int}[])]
        for i = 1:N-1
            append!(basis, [tuple([3i-2;3j+1], Vector{Int}[]) for j = i:N-1])
            # append!(basis, [tuple([3i-1;3j+2], Vector{Int}[]) for j = i:N-1])
            # append!(basis, [tuple([3i;3j+3], Vector{Int}[]) for j = i:N-1])
        end
        if d > 2
            push!(basis, tuple([1;5;9], Vector{Int}[]), tuple([1;6;8], Vector{Int}[]), tuple([2;4;9], Vector{Int}[]), tuple([2;6;7], Vector{Int}[]), tuple([3;4;8], Vector{Int}[]), tuple([3;5;7], Vector{Int}[]))
            push!(basis, tuple([1;11;15], Vector{Int}[]), tuple([1;12;14], Vector{Int}[]), tuple([2;10;15], Vector{Int}[]), tuple([2;12;13], Vector{Int}[]), tuple([3;10;14], Vector{Int}[]), tuple([3;11;13], Vector{Int}[]))
        end
    else
        basis = [tuple([3i-2], Vector{Int}[]) for i = 1:N]
        if d > 2
            push!(basis, tuple([1;4;7], Vector{Int}[]), tuple([1;5;8], Vector{Int}[]), tuple([2;4;8], Vector{Int}[]), tuple([2;5;7], Vector{Int}[]), tuple([1;6;9], Vector{Int}[]), tuple([3;4;9], Vector{Int}[]), tuple([3;6;7], Vector{Int}[]))
            push!(basis, tuple([1;10;13], Vector{Int}[]), tuple([1;11;14], Vector{Int}[]), tuple([2;10;14], Vector{Int}[]), tuple([2;11;13], Vector{Int}[]), tuple([1;12;15], Vector{Int}[]), tuple([3;10;15], Vector{Int}[]), tuple([3;12;13], Vector{Int}[]))
        end
        for a in edges
            append!(basis, [tuple([3*a[1]-1;3*a[2]], Vector{Int}[]), tuple([3*a[1];3*a[2]-1], Vector{Int}[])])
        end
        for a in triples
            append!(basis, [tuple([3*a[1]-1;3*a[2]], Vector{Int}[]), tuple([3*a[1]-1;3*a[3]], Vector{Int}[]), tuple([3*a[2]-1;3*a[3]], Vector{Int}[]), tuple([3*a[1];3*a[2]-1], Vector{Int}[]), tuple([3*a[1];3*a[3]-1], Vector{Int}[]), tuple([3*a[2];3*a[3]-1], Vector{Int}[])])
        end
    end
    return basis
end

function get_kagome_sbasis(N)
    basis = [[Vector{Int}[]]; [[[3i-2]] for i = 1:N]]
    return basis
end

function get_kagome_wbasis(N, d; label=1)
    if label == 1
        basis = []
    else
        if N == 5
            basis = [[1]]
        elseif N == 13
            basis = [[3i-2] for i = 1:5]
        elseif N == 27
            basis = [[3i-2] for i = 1:13]
        elseif N == 45
            basis = [[3i-2] for i = 1:27]
        else
            @error "Wrong number of sites!"
        end
    end
    return basis
end
