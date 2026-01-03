using SpectralGap

# 1D Ising model with transversal field
N = 18
H = ncpoly([[3*[i;i+1] for i = 1:N-1]; [[3*i-2] for i = 1:N]], [-ones(N-1); 0.5*ones(N)])
gamma = 0.29
d = 3
certify_gap(N, H, gamma, d, QUIET=false)
