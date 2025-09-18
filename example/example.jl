using SpectralGap

# 1D Ising model with transversal field
N = 3 # number of sites
H = ncpoly([[[i;i+1] for i = 1:N-1]; [[1;N]]], -ones(N)) # define the Hamiltonian
gamma = 0.5
d = 2
@time status = certify_gap(N, H, gamma, d, QUIET=false)

N = 2
H = ncpoly([[3;6]], [-1])
gamma = 0.5
d = 2
@time status = certify_gap(N, H, gamma, d, QUIET=false)


N = 2
H = ncpoly([[3;6], [1], [4]], [-1, 2, 2])
gamma = 1.3
d = 2
@time status = certify_gap(N, H, gamma, d, QUIET=false)

N = 3
H = ncpoly([[3;6], [6;9], [3;9], [1], [4], [7]], [-1, -1, -1, 2, 2, 2])
gamma = 1.5
d = 2
@time status = certify_gap(N, H, gamma, d, QUIET=false)

N = 4
H = ncpoly([[3;6], [6;9], [9;12], [3;12], [1], [4], [7], [10]], [-1, -1, -1, -1, 2, 2, 2, 2])
gamma = 0.5
d = 2
@time status = certify_gap(N, H, gamma, d, QUIET=false)

using LinearAlgebra

Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
N = 3
B = zeros(Int, 2^N, 2^N)
for i = 1:N-1
    ind = ones(Int, N)
    ind[i] = ind[i+1] = 4
    B -= real(kron(Pauli[ind]...))
end
ind = ones(Int, N)
ind[1] = ind[N] = 4
B -= real(kron(Pauli[ind]...))
for i = 1:N
    ind = ones(Int, N)
    ind[i] = 2
    B += 2*real(kron(Pauli[ind]...))
end
v = eigvals(B)
println(v[1:2])
println(v[2]-v[1])