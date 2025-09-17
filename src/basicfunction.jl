function get_basis(L, label, d; lattice="chain", extra=0, three_type=[1;1])
    basis = Vector{UInt16}[]
    
    return basis
end

# binary search in a sorted sequence
function bfind(A, l, a)
    low = 1
    high = l
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

# reduction to the normal form
function reduce1!(a::Vector{UInt16})
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
function reduce2!(a::Vector{UInt16}; realify=false)
    la = length(a)
    flag = 1
    coef = 1
    while flag == 1
        ind = findfirst(x -> a[x] != a[x+1] && ceil(Int, a[x]/3) == ceil(Int, a[x+1]/3), 1:la-1)
        if ind !== nothing
            s = mod.(a[ind:ind+1], 3)
            if s == [1, 2]
                a[ind] += UInt16(2)
                coef *= im
            elseif s == [0, 2]
                a[ind] -= UInt16(2)
                coef *= -im
            elseif s == [2, 1] || s == [1, 0]
                a[ind] += UInt16(1)
                coef *= -im
            else
                a[ind] -= UInt16(1)
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
function reduce3!(a::Vector{UInt16})
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
function isz(a::Vector{UInt16})
    return any(i->isodd(count(isequal(i), mod.(a,3))), 0:2)
end

# reduction w.r.t symmetries
function reduce4(a::Vector{UInt16}, L)
    l = length(a)
    if l > 0
        pa = Vector{UInt16}[]
        for i = 1:l
            ta = [a[i:end]; a[1:i-1] .+ 3*L] .- 3*(ceil(UInt16, a[i]/3) - 1)
            append!(pa, perm(ta))
            rta = reverse(ta)
            ma = 3*(ceil(UInt16, ta[end]/3) .- ceil.(UInt16, rta/3)) + smod.(rta, 3)
            append!(pa, perm(ma))
        end
        return findmin(pa)[1]
    else
        return a
    end
end

# implement all reductions
function reduce!(a::Vector{UInt16}; L=0, realify=false)
    reduce1!(a)
    reduce3!(a)
    a,coef = reduce2!(a, realify=realify)
    reduce3!(a)
    if isz(a)
        coef = 0
    else
        a = reduce4(a, L)
    end
    return a,coef
end

function rot(label)
    if label == 1
        return 2,3
    elseif label == 2
        return 3,1
    else
        return 1,2
    end
end

function smod(i, s)
    r = mod(i, s)
    return r == 0 ? s : r
end
