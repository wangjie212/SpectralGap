using SpectralGap

# 1D Ising model with transversal field
N = 3
H = ncpoly([[3*[i;i+1] for i = 1:N-1]; [[3*i-2] for i = 1:N]], [-ones(N-1); 2*ones(N)])
gamma = 1.5
d = 2
status,mmat = certify_gap(N, H, gamma, d, QUIET=false)
