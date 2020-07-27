include("intro.jl")

using DrWatson
include(srcdir("heightmaps.jl"))
include(srcdir("osn.jl"))
include(srcdir("exp_noise.jl"))


dims = (256, 256)
scale = 0.02
name = "Testeria"

map = init_heightmap(name, dims)

# cubic_noise_map!(map, scale, 8, 3.)
amortised_noise_map!(map)

island_modif = island_gradient(dims, 6. / 256, (128, 128), "quadratic", 2.0)

# apply_multiplier!(map, 0.1)
#
# apply_power!(map, 2.)
#
# apply_multiplier!(map, 30.)

apply_submap!(map, island_modif)

norm_map!(map)

apply_cutoff!(map, 0.3)

plot(1:256, 1:256, map.map, st = :heatmap, c = :greys)

# sns = surface_normals(map)
# plot(1:256, 1:256, sns[:,:,1], st = :heatmap, c = :greys)
