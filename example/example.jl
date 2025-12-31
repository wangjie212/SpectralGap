using SpectralGap

# 1D Ising model with transversal field
N = 5
H = ncpoly([[3*[i;i+1] for i = 1:N-1]; [[3*i-2] for i = 1:N]], [-ones(N-1); 2*ones(N)])
gamma = 2.87
d = 3
status,basis,mmat,tsupp = certify_gap(N, H, gamma, d, QUIET=false)


findfirst(a -> a == tuple(Int[], [[3*i-2], [3*j+1]]), tsupp)

any(a -> a == [[1], [4;10]], tsupp)

using Graphs

k = 2
G = SimpleGraph(size(mmat[k],1))
for i = 1:size(mmat[k],1), j = i+1:size(mmat[k],1)
    if abs(mmat[k][i,j]) > 1e-6
        add_edge!(G, i, j)
    end
end
blocks = connected_components(G)
println(length.(blocks))
