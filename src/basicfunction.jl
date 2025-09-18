mutable struct ncpoly
    supp::Vector{Vector{Int}}
    coe::Vector{Float64}
end

function get_wbasis(N, d)
    basis = [Int[]]
    for i = 1:3N
        push!(basis, [i])
    end
    if d > 1
        for i = 1:N-1, j = i+1:N, s = 1:3, t = 1:3
            push!(basis, [3*(i-1)+s,3*(j-1)+t])
        end
    end
    if d > 2
        for i = 1:N-2, j = i+1:N-1, k = j+1:N, s = 1:3, t = 1:3, u = 1:3
            push!(basis, [3*(i-1)+s,3*(j-1)+t,3*(k-1)+u])
        end
    end
    if d > 3
        for i = 1:N-3, j = i+1:N-2, k = j+1:N-1, l = k+1:N, s = 1:3, t = 1:3, u = 1:3, v = 1:3
            push!(basis, [3*(i-1)+s,3*(j-1)+t,3*(k-1)+u,3*(l-1)+v])
        end
    end
    return basis
end

function get_sbasis(wbasis, d)
    sbasis = [Int[]]
    i = 0
    temp = Int[]
    while i < d+1
        if sum(temp) == length(wbasis)*i
            temp = ones(Int, i+1)
            if i < d && sum(length.(wbasis[temp])) <= d
                push!(sbasis, temp)
            end
            i += 1
        else
            temp2 = copy(temp)
            j = temp[1]
            ind = findfirst(x->temp[x]!=j, 1:length(temp))
            if ind === nothing
                ind = length(temp)+1
            end
            if j != 1
                temp2[1:ind-2] = ones(Int, ind-2)
            end
            temp2[ind-1] = j+1
            temp = temp2
            if sum(length.(wbasis[temp])) <= d
                push!(sbasis, temp)
            end
        end
    end
    return sbasis
end

function get_tbasis(N, d, ptsupp)
    wbasis = get_wbasis(N, d)
    ind = findfirst(i->length(ptsupp[i])>d, 1:length(ptsupp))
    inx = ind !== nothing ? ind - 1 : length(ptsupp)
    sbasis = get_sbasis(ptsupp[1:inx], d)
    tbasis = Vector{Int}[]
    for i = 1:length(wbasis), j = 1:length(sbasis)
        if length(wbasis[i]) + sum(length.(ptsupp[sbasis[j]])) <= d
            push!(tbasis, [i,j])
        end
    end
    return tbasis,wbasis,sbasis
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
    if realify == true && !isreal(coef)
        coef = imag(coef)
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
    return any(i->isodd(count(isequal(i), mod.(a,3))), 0:2)
end

# reduction w.r.t symmetries
function reduce4(a::Vector{Int}, L)
    l = length(a)
    if l > 0
        pa = Vector{Int}[]
        for i = 1:l
            ta = [a[i:end]; a[1:i-1] .+ 3*L] .- 3*(ceil(Int, a[i]/3) - 1)
            append!(pa, perm(ta))
            rta = reverse(ta)
            ma = 3*(ceil(Int, ta[end]/3) .- ceil.(Int, rta/3)) + smod.(rta, 3)
            append!(pa, perm(ma))
        end
        return findmin(pa)[1]
    else
        return a
    end
end

# implement all reductions
function reduce!(a::Vector{Int}; realify=false)
    reduce1!(a)
    reduce3!(a)
    a,coef = reduce2!(a, realify=realify)
    reduce3!(a)
    # if isz(a)
    #     coef = 0
    # else
    #     a = reduce4(a, L)
    # end
    return a,coef
end

function smod(i, s)
    r = mod(i, s)
    return r == 0 ? s : r
end

function state_reduce(word1, word2, ptsupp)
    if isempty(word1)
        ind = Int[]
    else
        ind = bfind(ptsupp, length(ptsupp), word1, lt=isless_td)
    end
    return sort([word2; ind])
end

function isless_td(a, b)
    if length(a) < length(b)
        return true
    elseif length(a) > length(b)
        return false
    else
        return a < b
    end
end

# binary search in a sorted sequence
function bfind(A, l, a; lt=isless, rev=false)
    low = 1
    high = l
    while low <= high
        mid = Int(ceil(1/2*(low+high)))
        if isequal(A[mid], a)
           return mid
        elseif lt(A[mid], a)
            if rev == false
                low = mid+1
            else
                high = mid-1
            end
        else
            if rev == false
                high = mid-1
            else
                low = mid+1
            end
        end
    end
    return nothing
end
