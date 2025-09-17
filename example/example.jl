using SpectralGap

# 1D Ising model with transversal field
supp = Vector{UInt16}[[3;6]]
coe = [-1]
gamma = 1/2
d = 2
@time status = SpectralGap(supp, coe, gamma, d, QUIET=false)


supp = Vector{UInt16}[[3;6], [1], [4]]
coe = [-1, 2, 2]
gamma = 1/2
d = 2
@time status = SpectralGap(supp, coe, gamma, d, QUIET=false)