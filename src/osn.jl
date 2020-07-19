function generate_noisemap(dims::NTuple{2,Int64}, scale::Float64 = 0.2, octaves::Int64 = 5, falloff::Float64 = 2.0, SEED::Int64 = 1500493948)
   for i in dims
      if i < 1
         DomainError("Map size not positive")
      end
   end

   RND_A = 134775813
   RND_B = 1103515245

   #const SEED = 1500493948

   MAX_VAL = 2^31 - 1
   periodx = MAX_VAL
   periody = MAX_VAL

   MASK = 2^32 - 1

   function interpolate(a, b, c, d, x)
      p = (d - c) - (a - b)
      return (x * (x * (x * p + ((a - b) - p)) + (c - a)) + b)
   end

   function randomize(x, y, seed = SEED)
      return (((((x ⊻ y) * RND_A) ⊻ (seed + x)) * (((RND_B * x) << 16) ⊻ (RND_B * y) - RND_A)) & MASK) / MASK
   end

   function tile(coordinate, period = MAX_VAL)
      return coordinate % period
   end

   function sample1d(x, seed = SEED, octave = 1)
      xi = floor(Int64, x / octave)
      lerp = x / octave - xi
      return interpolate(
         randomize(tile(xi - 1, periodx), 0, seed),
         randomize(tile(xi, periodx), 0, seed),
         randomize(tile(xi + 1, periodx), 0, seed),
         randomize(tile(xi + 2, periodx), 0, seed),
         lerp,
      ) * 0.5 + 0.25
   end

   function sample2d(x, y, seed = SEED, octave = 1)
      xi::Int64 = floor(Int64, x / octave)
      lerpx::Float64 = x / octave - xi
      yi::Int64 = floor(Int64, y / octave)
      lerpy::Float64 = y / octave - yi

      xSamples = [0.0 0.0 0.0 0.0]

      for iin = 1:4
         ii = iin
         xSamples[iin] = interpolate(
            randomize(tile(xi - 1, periodx), tile(yi - 1 + ii, periody), seed),
            randomize(tile(xi, periodx), tile(yi - 1 + ii, periody), seed),
            randomize(tile(xi + 1, periodx), tile(yi - 1 + ii, periody), seed),
            randomize(tile(xi + 2, periodx), tile(yi - 1 + ii, periody), seed),
            lerpx,
         )
      end

      return interpolate(xSamples..., lerpy) * 0.5 + 0.25
   end

   tarr = zeros(dims[1], dims[2])
   dims = size(tarr)

   falloff_norm = 1.0 / sum(falloff .^ -(1:octaves))

   for oct = 1:octaves
      for col = 1:dims[2], row = 1:dims[1]
         tarr[row, col] += falloff_norm * sample2d(row * scale, col * scale, SEED, 1) * falloff^(-oct)
      end
   end

   return tarr
end

# dimm = 1000
# pn = generate_noisemap((dimm, dimm), 0.02)
#
# plot(1:dimm, 1:dimm, pn, st = :heatmap, c = :grays)
