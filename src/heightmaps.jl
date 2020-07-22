using Random
include("osn.jl")
include("exp_noise.jl")

mutable struct heightmap
    name::String
    seed::Int64
    map::Array{Float64,2}
end

function init_heightmap(name::String, dims::NTuple{2,Int64})
    for dim in dims
        @assert dim > 0
    end
    seed = rand(Int64)
    map = zeros(Float64, dims[1], dims[2])
    rv = heightmap(name, seed, map)
    return rv
end


function island_gradient(dims::NTuple{2,Int64}, sharpness::Float64, origin::NTuple{2,Int64}, method::String = "exponential", power::Float64 = 2.0)
    for dim in dims
        @assert dim > 0
    end
    grad_map = zeros(Float64, dims[1], dims[2])
    if method == "exponential"
        for col = 1:dims[2], row = 1:dims[1]
            grad_map[row, col] = exp(-(abs(sharpness * (col - origin[2]))^power + abs(sharpness * (row - origin[1]))^power))
        end
    elseif method == "quadratic"
        for col = 1:dims[2], row = 1:dims[1]
            grad_map[row, col] = - sqrt((col - origin[2])^2 + (row - origin[1])^2)
        end
    elseif method == "square"
        for col = 1:dims[2], row = 1:dims[1]
            grad_map[row, col] = max((col - origin[2])^2, (row - origin[1])^2)
        end
    else
        ArgumentError("No such method exists")
    end
    norm_map!(grad_map)

    return grad_map
end

island_gradient(dims::NTuple{2,Int64}, origin::NTuple{2,Int64}) = island_gradient(dims::NTuple{2,Int64}, 5.5 / maximum(dims), origin::NTuple{2,Int64})

function n_island_gradient(dims::NTuple{2,Int64}, n_islands::Int64) end

function cubic_noise_map!(hm::heightmap, scale::Float64 = 0.2, octaves::Int64 = 5, falloff::Float64 = 2.0)
    dims = size(hm.map)
    hm.map = generate_noisemap(dims, scale, octaves, falloff, hm.seed)
    return hm
end

function xorshift(x::Int64, s1::Int64, s2::Int64, s3::Int64)::Int64
    o = x
    o ⊻= (o << s1)
    o ⊻= (o >> s2)
    o ⊻= (o << s3)
    return o
end

function exp_noise_map!(hm::heightmap, octave_0::Int64 = 1, octave_1::Int64 = 4, falloff::Float64 = 2.0, omega::Float64 = 0.5, coords::NTuple{2,Int64} = (1, 1))
    @assert 0.0 <= omega <= 1.0
    @assert 0 < octave_0 <= octave_1
    @assert falloff > 0.0
    dims = size(hm.map)
    @assert dims[1] == dims[2]
    @assert floor(log2(dims[1])) == log2(dims[1])
    seed_1 = xorshift(hm.seed, 8, 31, 17)
    seed_2 = xorshift(hm.seed, 13, 17, 5)
    hm.map = AMN_EXP_2D_grid(coords, seed_1, seed_2, omega, 1.0 / falloff, octave_0, octave_1, dims[1])
    return hm
end

function amortised_noise_map!(hm::heightmap, octave_0::Int64 = 3, octave_1::Int64 = 6, falloff::Float64 = 2.0, coords::NTuple{2,Int64} = (1, 1))
    @assert 0 < octave_0 <= octave_1
    @assert falloff > 0.0
    dims = size(hm.map)
    @assert dims[1] == dims[2]
    @assert floor(log2(dims[1])) == log2(dims[1])
    seed_1 = xorshift(hm.seed, 8, 31, 17)
    seed_2 = xorshift(hm.seed, 13, 17, 5)
    hm.map = AMN_EXP_2D_grid(coords, seed_1, seed_2, 0.5, 1.0 / falloff, octave_0, octave_1, dims[1])
    return hm
end

function apply_modifier!(hm::heightmap, mod::Array{Float64,2})
    @assert size(hm.map) == size(mod)
    hm.map .*= mod
    return hm
end

function apply_multiplier!(hm::heightmap, mult::Float64)
    hm.map *= mult
    return hm
end

function apply_submap!(hm::heightmap, sm::Array{Float64,2})
    @assert size(hm.map) == size(sm)
    hm.map = hm.map - (ones(size(sm)) - sm)
    return hm
end

function apply_power!(hm::heightmap, power::Float64)
    hm.map .^= power
    return hm
end

function apply_cutoff!(hm::heightmap, cutoff::Float64)
    hm.map = max.(hm.map, cutoff)
    return hm
end

function norm_map!(hm::heightmap)
    h_min = minimum(hm.map)
    hm.map .-= h_min
    h_max = maximum(hm.map)
    hm.map ./= h_max
    return hm
end

function norm_map!(hm::Array{Float64,2})
    h_min = minimum(hm)
    hm .-= h_min
    h_max = maximum(hm)
    hm ./= h_max
    return hm
end
