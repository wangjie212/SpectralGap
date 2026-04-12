module SpectralGap

using MathOptInterface
using JuMP
using MosekTools
using LinearAlgebra

export certify_Ising_gap, certify_Heisenberg_kagome_gap, ncpoly

include("basicfunction.jl")
include("sdp.jl")

end