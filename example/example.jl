using SpectralGap

# 1D Ising model with transversal field
N = 15
g = 0.5
H = ncpoly( [[3*[i;i+1] for i =1:N-1]; [[3i-2] for i = 1:N]], [-ones(N-1); g*ones(N)])
d = 3
ub = 0.15
lb = 0.15
while ub - lb > 1e-2
    gamma = (ub + lb)/2
    flag = certify_Ising_gap(N, H, gamma, d, QUIET=false)
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


N = 9
g = 0.5
H = ncpoly( [[3*[i;i+1] for i =1:N-1]; [[3i-2] for i = 1:N]], [-ones(N-1); g*ones(N)])
d = 3
ub = 2.90
lb = 2.80
while ub - lb > 1e-2
    gamma = (ub + lb)/2
    flag = certify_Ising_gap_nosignsymmetry(N, H, gamma, d, QUIET=false)
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

a = [0.5, 1.5, 2.5, 3.5]
b = [4.00, 2.52, 2.29, 2.19]
c = [3.00, 1.66, 1.37, 1.26]
d = [2.00, 0.93, 0.60, 0.47]
e = [1.00, 0.39, 0.19, 0.10]
p = plot(a, b, dpi=600, shape=:circle, legendfontsize=10, label=["g = 2.0"], xlabel="relaxation level (L, d)", xlims = (0.4,3.6), xdiscrete_values=[(2,1),(3,2),(4,3),(5,4)], ylabel="upper bounds on the spectral gap")
plot!(p, a, c, dpi=600, shape=:star, label=["g = 1.5"])
plot!(p, a, d, dpi=600, shape=:dtriangle, label=["g = 1.0"])
plot!(p, a, e, dpi=600, shape=:hexagon, label=["g = 0.5"])
plot!(p, range(0.4, 3.5, length=100), 2*ones(100), line = (:steppre, :dot, false, 0.5, 1, :red), label=nothing)
plot!(p, range(0.4, 3.5, length=100), ones(100), line = (:steppre, :dot, false, 0.5, 1, :red), label=nothing)
plot!(p, range(0.4, 3.5, length=100), zeros(100), line = (:steppre, :dot, false, 0.5, 1, :red), label=nothing)

savefig("d:/Programs/SpectralGap/Ising_spectral_gap.png")


# Kagome lattice Heisenberg model
N = 5
triples = [[1,2,3], [1,4,5]]
edges = []
N = 13
triples0 = [[1,2,3], [1,4,5], [2,6,7], [3,8,9], [4,10,11], [5,12,13]]
edges0 = []
N = 27
triples = [[1,2,3], [1,4,5], [2,6,7], [3,8,9], [4,10,11], [5,12,13], [6,14,27], [7,15,16], [8,17,18], [9,19,20], [10,20,21], [11,22,23], [12,24,25], [13,26,27]]
edges = [[16, 17], [23, 24]]
N = 45
triples = [[1,2,3], [1,4,5], [2,6,7], [3,8,9], [4,10,11], [5,12,13], [6,14,27], [7,15,16], [8,17,18], [9,19,20], [10,20,21], [11,22,23], [12,24,25], [13,26,27],
[14, 28, 29], [15, 30, 31], [16, 17, 32], [18, 33, 34], [19, 35, 36], [21, 37, 38], [22, 39, 40], [23, 24, 41], [25, 42, 43], [26, 44, 45]]
edges = [[29, 30], [34, 35], [38, 39], [43, 44]]

H = ncpoly(vcat([[[3*a[1]-2;3*a[2]-2], [3*a[1]-1;3*a[2]-1], [3*a[1];3*a[2]], [3*a[1]-2;3*a[3]-2], [3*a[1]-1;3*a[3]-1], [3*a[1];3*a[3]], [3*a[2]-2;3*a[3]-2], [3*a[2]-1;3*a[3]-1], [3*a[2];3*a[3]]] for a in triples]...), ones(9*length(triples)))
d = 3
ub = 1.21
lb = 1.21
while ub - lb > 1e-2
    gamma = (ub + lb)/2
    flag = certify_Heisenberg_kagome_gap(N, H, triples, edges, triples0, edges0, gamma, d, QUIET=false)
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


H = ncpoly(vcat([[[3*a[1]-2;3*a[2]-2], [3*a[1]-1;3*a[2]-1], [3*a[1];3*a[2]], [3*a[1]-2;3*a[3]-2], [3*a[1]-1;3*a[3]-1], [3*a[1];3*a[3]], [3*a[2]-2;3*a[3]-2], [3*a[2]-1;3*a[3]-1], [3*a[2];3*a[3]]] for a in triples]...), ones(9*length(triples)))
obj = ncpoly([[1;5;9], [2;6;7], [3;4;8], [1;6;8], [3;5;7], [2;4;9]], 1/8*[1; 1; 1; -1; -1; -1])
d = 2
ub = 1.31
lb = 1.31
while ub - lb > 1e-2
    gamma = (ub + lb)/2
    flag,v = certify_Heisenberg_kagome_gap_nosignsymmetry(N, H, triples, edges, gamma, d, QUIET=false)
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
c = [2.01, 1.28, 1.15]
d = [2.00, 1.31]
p = plot(a, b, dpi=600, shape=:circle, legendfontsize=10, label=["all symmetries, d = 2"], xlabel="L", ylabel="upper bounds on the spectral gap")
plot!(p, a[1:3], c, dpi=600, shape=:star, label=["all symmetries, d = 3"])
plot!(p, a[1:2], d, dpi=600, shape=:dtriangle, label=["pi-spin rotation symmetry, d = 2"])

savefig("d:/Programs/SpectralGap/Kagome.png")