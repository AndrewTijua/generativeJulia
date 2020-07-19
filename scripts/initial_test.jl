using DrWatson
include(srcdir("heightmaps.jl"))
include(srcdir("osn.jl"))

using Plots
plotlyjs()

dims = (256, 256)
scale = 0.02
name = "Testeria"

map = init_heightmap(name, dims)

cubic_noise_map!(map, scale, 8, 3.)

island_modif = island_gradient(dims, 5.5 / 256, (128, 128), 3.3)

# apply_multiplier!(map, 0.1)
#
# apply_power!(map, 2.)
#
# apply_multiplier!(map, 30.)

apply_submap!(map, island_modif)

norm_map!(map)

# apply_cutoff!(map, 0.2)

plot(1:256, 1:256, map.map, st = :surface, c = :greys)
