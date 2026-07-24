function posepsd6!(model, cons, tsupp, L; lattice="Ising", identify_zeros=true, symmetry=true)
    Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    blocks = [Vector(1:2^6)]
    pos = [@variable(model, [1:length(block), 1:length(block)], PSD) for block in blocks]
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3
        ind = [i,j,k,l,s,t]
        inx = ind .!= 0
        bi = 3*(Vector(1:6)[inx] .- 1) + ind[inx]
        if !identify_zeros || !isz(bi, model=lattice)
            if symmetry == true
                if lattice == "Ising"
                    bi = reduce_mirror(bi, L)
                else
                    bi = reduce_perm(bi)
                end
            end
            Locb = bfind(tsupp, [bi])
            tp = real(kron(Pauli[ind.+1]...))
            for (p, block) in enumerate(blocks)
                @inbounds add_to_expression!(cons[Locb], sum(tp[block, block].*pos[p]))
            end
        end
    end
end

function posepsd7!(model, cons, tsupp, L; lattice="Ising", identify_zeros=true, symmetry=true)
    Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    blocks = [Vector(1:2^7)]
    pos = [@variable(model, [1:length(block), 1:length(block)], PSD) for block in blocks]
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3
        ind = [i,j,k,l,s,t,u]
        inx = ind .!= 0
        bi = 3*(Vector(1:7)[inx] .- 1) + ind[inx]
        if !identify_zeros || !isz(bi, model=lattice)
            if symmetry == true
                if lattice == "Ising"
                    bi = reduce_mirror(bi, L)
                else
                    bi = reduce_perm(bi)
                end
            end
            Locb = bfind(tsupp, [bi])
            tp = real(kron(Pauli[ind.+1]...))
            for (p, block) in enumerate(blocks)
                @inbounds add_to_expression!(cons[Locb], sum(tp[block, block].*pos[p]))
            end
        end
    end
end

function posepsd8!(model, cons, tsupp, L; lattice="Ising", identify_zeros=true, symmetry=true)
    Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    blocks = [Vector(1:2^8)]
    pos = [@variable(model, [1:length(block), 1:length(block)], PSD) for block in blocks]
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3
        ind = [i,j,k,l,s,t,u,v]
        inx = ind .!= 0
        bi = 3*(Vector(1:8)[inx] .- 1) + ind[inx]
        if !identify_zeros || !isz(bi, model=lattice)
            if symmetry == true
                if lattice == "Ising"
                    bi = reduce_mirror(bi, L)
                else
                    bi = reduce_perm(bi)
                end
            end
            Locb = bfind(tsupp, [bi])
            tp = real(kron(Pauli[ind.+1]...))
            for (p, block) in enumerate(blocks)
                @inbounds add_to_expression!(cons[Locb], sum(tp[block, block].*pos[p]))
            end
        end
    end
end

function posepsd9!(model, cons, tsupp)
    Pauli = Matrix{Complex{Int8}}[[1 0; 0 1], [0 1; 1 0], [0 -im; im 0], [1 0; 0 -1]]
    blocks = [[1], [2, 3, 5, 9, 17, 33, 65, 129, 257], [4, 6, 7, 10, 11, 13, 18, 19, 21, 25, 34, 35, 37, 41, 49, 66, 67, 69, 73, 81, 97, 130, 131, 133, 137, 145, 161, 193, 258, 259, 261, 265, 273, 289, 321, 385], [8, 12, 14, 15, 20, 22, 23, 26, 27, 29, 36, 38, 39, 42, 43, 45, 50, 51, 53, 57, 68, 70, 71, 74, 75, 77, 82, 83, 85, 89, 98, 99, 101, 105, 113, 132, 134, 135, 138, 139, 141, 146, 147, 149, 153, 162, 163, 165, 169, 177, 194, 195, 197, 201, 209, 225, 260, 262, 263, 266, 267, 269, 274, 275, 277, 281, 290, 291, 293, 297, 305, 322, 323, 325, 329, 337, 353, 386, 387, 389, 393, 401, 417, 449], [16, 24, 28, 30, 31, 40, 44, 46, 47, 52, 54, 55, 58, 59, 61, 72, 76, 78, 79, 84, 86, 87, 90, 91, 93, 100, 102, 103, 106, 107, 109, 114, 115, 117, 121, 136, 140, 142, 143, 148, 150, 151, 154, 155, 157, 164, 166, 167, 170, 171, 173, 178, 179, 181, 185, 196, 198, 199, 202, 203, 205, 210, 211, 213, 217, 226, 227, 229, 233, 241, 264, 268, 270, 271, 276, 278, 279, 282, 283, 285, 292, 294, 295, 298, 299, 301, 306, 307, 309, 313, 324, 326, 327, 330, 331, 333, 338, 339, 341, 345, 354, 355, 357, 361, 369, 388, 390, 391, 394, 395, 397, 402, 403, 405, 409, 418, 419, 421, 425, 433, 450, 451, 453, 457, 465, 481]]
    pos = [@variable(model, [1:length(block), 1:length(block)], PSD) for block in blocks]
    for i = 0:3, j = 0:3, k = 0:3, l = 0:3, s = 0:3, t = 0:3, u = 0:3, v = 0:3, w = 0:3
        ind = [i,j,k,l,s,t,u,v,w]
        if all(x->iseven(sum(ind .== x)), 1:3)
            inx = ind .!= 0
            bi = 3*(Vector(1:9)[inx] .- 1) + ind[inx]
            Locb = bfind(tsupp, [reduce_perm(bi)])
            tp = real(kron(Pauli[ind.+1]...))
            for (p, block) in enumerate(blocks)
                @inbounds add_to_expression!(cons[Locb], sum(tp[block, block].*pos[p]))
            end
        end
    end
end
