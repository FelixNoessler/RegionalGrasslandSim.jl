function get_allneighbours(; p)
    xs = [[x for x in 1:p.patch_xdim, _ in 1:p.patch_ydim]...]
    ys = [[y for _ in 1:p.patch_xdim, y in 1:p.patch_ydim]...]

    neighbours = Vector{Int64}[]
    for pa in Base.OneTo(p.npatches)
        x = xs[pa]
        y = ys[pa]
        push!(neighbours, get_neighbours(; x, y, p))
    end

    return neighbours
end

function get_neighbours(; x, y, p)
    neighbours = Int64[]
    if x - 1 > 0
        push!(neighbours, get_patchindex(x-1, y; p))
    end
    if x + 1 <= p.patch_xdim
        push!(neighbours, get_patchindex(x+1, y; p))
    end
    if y + 1 <= p.patch_ydim
        push!(neighbours, get_patchindex(x, y+1; p))
    end
    if y - 1 > 0
        push!(neighbours, get_patchindex(x, y-1; p))
    end

    return neighbours
end

function get_patchindex(x, y; p)
    return x + (y-1) * p.patch_xdim
end
