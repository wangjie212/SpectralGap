using SpectralGap

# 1D Ising model with transversal field
N = 5
H = ncpoly([[3*[i;i+1] for i = 1:N-1]; [[3*i-2] for i = 1:N]], [-ones(N-1); 2*ones(N)])
d = 3
ub = 3
lb = 2
while ub - lb > 1e-2
    gamma = (ub + lb)/2
    flag = certify_gap(N, H, gamma, d, QUIET=true)
    if flag == 1
        lb = gamma
    else
        ub = gamma
    end
    println([lb, ub])
end
flag = certify_gap(N, H, trunc(ub, digits=2), d, QUIET=true)
if flag == 1
    println("Upper Bound = $(trunc(ub, digits=2)+0.01)")
else
    println("Upper Bound = $(trunc(ub, digits=2))")
end
