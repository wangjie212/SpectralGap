using SpectralGap

# 1D Ising model with transversal field
N = 3 # number of sites
H = ncpoly([3*[i;i+1] for i = 1:N-1], -ones(N-1)) # define the Hamiltonian
gamma = 3.5
d = 2
@time status = certify_gap(N, H, gamma, d, QUIET=false)

N = 5
H = ncpoly([[3*[i;i+1] for i = 1:N-1]; [[3*i-2] for i = 1:N]], [-ones(N-1); 1.5*ones(N)])
gamma = 3.5
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
# ind = ones(Int, N)
# ind[1] = ind[N] = 4
# B -= real(kron(Pauli[ind]...))
for i = 1:N
    ind = ones(Int, N)
    ind[i] = 2
    B += 2*real(kron(Pauli[ind]...))
end
v = eigvals(B)
println(v[1:2])
println(v[2]-v[1])