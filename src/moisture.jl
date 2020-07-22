include(scriptsdir("intro.jl"))

using DrWatson
include(srcdir("heightmaps.jl"))
include(srcdir("osn.jl"))
include(srcdir("exp_noise.jl"))

mutable struct biomemap
    name::String
    seed::Int64
    height::Array{Float64,2}
    moisture::Array{Float64,2}
end

function convert_h2b(hm::heightmap, mm::Array{Float64,2})
    @assert size(mm) == size(hm.map)
    bm = biomemap(hm.name, hm.seed, hm.map, mm)
    return bm
end

function generate_rendertypes(bm::biomemap)
    dims = size(bm.height)

    types = [
        "stde" "glad" "trsf" "trrf"
        "tede" "glad" "tedf" "terf"
        "tede" "shld" "taig" "taig"
        "scla" "bala" "tund" "snow"
    ]

    water_lev = minimum(bm.height)

    rendertypes = Array{String,2}(undef, dims[1], dims[2])
    for col = 1:dims[2], row = 1:dims[1]
        m_offset = ceil(Int64, bm.moisture[row, col] * 4.0)
        h_offset = ceil(Int64, bm.height[row, col] * 4.0)
        if bm.height[row, col] == water_lev
            rendertypes[row, col] = "wate"
        else
            rendertypes[row, col] = types[m_offset, h_offset]
        end
    end
    return rendertypes
end
