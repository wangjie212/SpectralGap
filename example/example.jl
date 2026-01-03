using SpectralGap

# 1D Ising model with transversal field
N = 18
H = ncpoly([[3*[i;i+1] for i = 1:N-1]; [[3*i-2] for i = 1:N]], [-ones(N-1); 0.5*ones(N)])
gamma = 0.29
d = 3
status,basis,mmat = certify_gap(N, H, gamma, d, QUIET=false)

using Graphs

k = 1
G = SimpleGraph(size(mmat[k],1))
for i = 1:size(mmat[k],1), j = i+1:size(mmat[k],1)
    if abs(mmat[k][i,j]) > 1e-6
        add_edge!(G, i, j)
    end
end
blocks = connected_components(G)
println(length.(blocks))
