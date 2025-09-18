module SpectralGap

using JuMP
using MosekTools
using LinearAlgebra

export certify_gap, ncpoly

include("basicfunction.jl")
include("sdp.jl")

end