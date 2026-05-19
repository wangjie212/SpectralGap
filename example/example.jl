using SpectralGap

# 1D Ising model with transversal field
N = 5
g = 0.5
H = ncpoly( [[3*[i;i+1] for i =1:N-1]; [[3i-2] for i = 1:N]], [-ones(N-1); g*ones(N)])
d = 4
ub = 1
lb = 0
while ub - lb > 1e-2
    gamma = (ub + lb)/2
    flag = certify_Ising_gap(N, H, gamma, d, QUIET=true)
    if flag == 1
        lb = gamma
    else
        ub = gamma
    end
    println([lb, ub])
end
flag = certify_Ising_gap(N, H, trunc(ub, digits=2), d, QUIET=true)
if flag == 1
    println("Upper Bound = $(trunc(ub, digits=2)+0.01)")
else
    println("Upper Bound = $(trunc(ub, digits=2))")
end

# Only mirror symmerty
N = 5
g = 0.5
H = ncpoly( [[3*[i;i+1] for i =1:N-1]; [[3i-2] for i = 1:N]], [-ones(N-1); g*ones(N)])
d = 3
ub = 3.5
lb = 3.4
while ub - lb > 1e-2
    gamma = (ub + lb)/2
    flag = certify_Ising_gap_nosignsymmetry(N, H, gamma, d, QUIET=true)
    if flag == 1
        lb = gamma
    else
        ub = gamma
    end
    println([lb, ub])
end
flag = certify_Ising_gap_nosignsymmetry(N, H, trunc(ub, digits=2), d, QUIET=true)
if flag == 1
    println("Upper Bound = $(trunc(ub, digits=2)+0.01)")
else
    println("Upper Bound = $(trunc(ub, digits=2))")
end


using Plots

a = Vector(3:2:19)
b = [4.01, 3.30, 2.88, 2.63, 2.48, 2.38, 2.30, 2.25, 2.22, 2.19, 2.17, 2.15, 2.14, 2.13, 2.12, 2.11, 2.11, 2.10]
c = [3.01, 2.36, 1.98, 1.74, 1.59, 1.49, 1.42, 1.37, 1.34, 1.31, 1.29, 1.28, 1.27, 1.26, 1.26, 1.25, 1.25, 1.25]
d = [2.01, 1.47, 1.14, 0.94, 0.83, 0.77, 0.73, 0.71, 0.69, 0.68, 0.68, 0.67, 0.67, 0.66, 0.66, 0.66, 0.66, 0.65]
e = [1.01, 0.65, 0.44, 0.36, 0.32, 0.30, 0.29, 0.29, 0.28, 0.28, 0.28, 0.28, 0.27, 0.27, 0.27, 0.27, 0.27, 0.27]
p = plot(a, b[a.-2], dpi=600, shape=:circle, label=["g = 2.0"], xlabel="L", xlims = (2,20), ylabel="upper bounds on the spectral gap")
plot!(p, a, c[a.-2], dpi=600, shape=:star, label=["g = 1.5"])
plot!(p, a, d[a.-2], dpi=600, shape=:dtriangle, label=["g = 1.0"])
plot!(p, a, e[a.-2], dpi=600, shape=:hexagon, label=["g = 0.5"])
plot!(p, range(2, 19, length=100), 2*ones(100), line = (:steppre, :dot, false, 0.5, 1, :red), label=nothing)
plot!(p, range(2, 19, length=100), ones(100), line = (:steppre, :dot, false, 0.5, 1, :red), label=nothing)
plot!(p, range(2, 19, length=100), zeros(100), line = (:steppre, :dot, false, 0.5, 1, :red), label=nothing)

savefig("d:/Programs/SpectralGap/Ising_spectral_gap.png")


# Kagome lattice Heisenberg model
N = 5
triples0 = [[1,2,3], [1,4,5]]
edges0 = []
N = 13
triples = [[1,2,3], [1,4,5], [2,6,7], [3,8,9], [4,10,11], [5,12,13]]
edges = []
N = 27
triples = [[1,2,3], [1,4,5], [2,6,7], [3,8,9], [4,10,11], [5,12,13], [6,14,27], [7,15,16], [8,17,18], [9,19,20], [10,20,21], [11,22,23], [12,24,25], [13,26,27]]
edges = [[16, 17], [23, 24]]
N = 45
triples = [[1,2,3], [1,4,5], [2,6,7], [3,8,9], [4,10,11], [5,12,13], [6,14,27], [7,15,16], [8,17,18], [9,19,20], [10,20,21], [11,22,23], [12,24,25], [13,26,27],
[14, 28, 29], [15, 30, 31], [16, 17, 32], [18, 33, 34], [19, 35, 36], [21, 37, 38], [22, 39, 40], [23, 24, 41], [25, 42, 43], [26, 44, 45]]
edges = [[29, 30], [34, 35], [38, 39], [43, 44]]

H = ncpoly(vcat([[[3*a[1]-2;3*a[2]-2], [3*a[1]-1;3*a[2]-1], [3*a[1];3*a[2]], [3*a[1]-2;3*a[3]-2], [3*a[1]-1;3*a[3]-1], [3*a[1];3*a[3]], [3*a[2]-2;3*a[3]-2], [3*a[2]-1;3*a[3]-1], [3*a[2];3*a[3]]] for a in triples]...), ones(9*length(triples)))
d = 3
ub = 2.1
lb = 2
while ub - lb > 1e-2
    gamma = (ub + lb)/2
    flag = certify_Heisenberg_kagome_gap(N, H, triples, edges, triples0, edges0, gamma, d, QUIET=true)
    if flag == 1
        lb = gamma
    else
        ub = gamma
    end
    println([lb, ub])
end
flag = certify_Heisenberg_kagome_gap(N, H, triples, edges, triples0, edges0, trunc(ub, digits=2), d, QUIET=true)
if flag == 1
    println("Upper Bound = $(trunc(ub, digits=2)+0.01)")
else
    println("Upper Bound = $(trunc(ub, digits=2))")
end


# Only pi-spin symmetry
N = 5
triples = [[1,2,3], [1,4,5]]
edges = []
N = 13
triples = [[1,2,3], [1,4,5], [2,6,7], [3,8,9], [4,10,11], [5,12,13]]
edges = []

H = ncpoly(vcat([[[3*a[1]-2;3*a[2]-2], [3*a[1]-1;3*a[2]-1], [3*a[1];3*a[2]], [3*a[1]-2;3*a[3]-2], [3*a[1]-1;3*a[3]-1], [3*a[1];3*a[3]], [3*a[2]-2;3*a[3]-2], [3*a[2]-1;3*a[3]-1], [3*a[2];3*a[3]]] for a in triples]...), ones(9*length(triples)))
obj = ncpoly([[1;5;9], [2;6;7], [3;4;8], [1;6;8], [3;5;7], [2;4;9]], 1/8*[1; 1; 1; -1; -1; -1])
d = 2
ub = 2.1
lb = 1.9
while ub - lb > 1e-2
    gamma = (ub + lb)/2
    flag,v = certify_Heisenberg_kagome_gap_nosignsymmetry(N, H, triples, edges, gamma, d, QUIET=true)
    if flag == 1
        lb = gamma
        println([lb, ub, v])
    else
        ub = gamma
    end
end
println("Upper Bound = $ub")


using Plots

a = Vector(1:4)
b = [2.01, 1.31, 1.24, 1.22]
p = plot(a, b, dpi=600, shape=:circle, label=["kagome"], xlabel="L", ylabel="upper bounds on the spectral gap")

savefig("d:/Programs/SpectralGap/Kagome_spectral_gap.png")