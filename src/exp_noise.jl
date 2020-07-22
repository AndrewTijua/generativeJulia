using Random
using DrWatson
# using Plots
# gr()

s_curve2(t) = (t * t * t * (10.0 + 3.0 * t * (2.0 * t - 5.0)))
s_curve(x) = x^2 * (3 - 2 * x)
lerp(t, a, b) = a + t * (b - a)

clip(x, a, b) = min(max(x, a), b)

rhash(a, b) = hash(hash(a) + hash(b))
rhash(a, b, seed) = hash(hash(a) + hash(b) + hash(seed))

function UniformHash(x::Int64, max::Int64)
    return (x + 1.0) / (max + 2.0)
end

function UniformHash(x::UInt64, max::Int64)
    return (x + 1.0) / (max + 2.0)
end

function UniformHash(x::UInt64, max::UInt64)
    return (x + 1.0) / (max + 2.0)
end

function ExpHash(x::Int64, max::Int64)
    scale::Float64 = 1. / log(0.5 * (max + 2.0))
    return (-scale * log(0.5 * x + 1.0) + 1.0)
end

function ExpHash(x::UInt64, max::Int64)
    scale::Float64 = 1. / log(0.5 * (max + 2.0))
    return (-scale * log(0.5 * x + 1.0) + 1.0)
end

function ExpHash(x::UInt64, max::UInt64)
    scale::Float64 = 1. / log(0.5 * (max + 2.0))
    return (-scale * log(0.5 * x + 1.0) + 1.0)
end

function ExpHash(x::Int64, y::Int64, m::Int64, omega::Float64)
    omega = clip(omega, 0.0, 1.0)

    return ((UniformHash(y, m) < omega) ? UniformHash(x, m) : ExpHash(x, m))
end

function ExpHash(x::UInt64, y::UInt64, m::UInt64, omega::Float64)
    omega = clip(omega, 0.0, 1.0)

    return ((UniformHash(y, m) < omega) ? UniformHash(x, m) : ExpHash(x, m))
end

function ExpHash(x::UInt64, y::UInt64, m::Int64, omega::Float64)
    omega = clip(omega, 0.0, 1.0)

    return ((UniformHash(y, m) < omega) ? UniformHash(x, m) : ExpHash(x, m))
end

function init_spline_table(n::Int64)
    s_table = zeros(n)
    for i = 1:n
        t = (i - 1) / n
        s_table[i] = s_curve2(t)
    end
    return s_table
end

function spline_table_inplace!(st, n::Int64)
    for i = 1:n
        t = (i - 1) / n
        st[i] = s_curve2(t)
    end
    return st
end

function FillDn!(t::Array{Float64,1}, s::Float64, n::Int64)
    d = -s / n
    t[n] = d
    for i = (n-1):-1:1
        t[i] = t[i+1] + d
    end
    return t
end

function FillUp!(t::Array{Float64,1}, s::Float64, n::Int64)
    d = s / n
    t[1] = 0.0
    for i = 2:n
        t[i] = t[i-1] + d
    end
    return t
end


function init_edge_tables_exp(x0::Int64, y0::Int64, n::Int64, omega::Float64, seed1, seed2, igt)
    b00 = rhash(x0, y0, seed1)
    b01 = rhash(x0, y0 + 1, seed1)
    b10 = rhash(x0 + 1, y0, seed1)
    b11 = rhash(x0 + 1, y0 + 1, seed1)

    max = typemax(UInt64)

    m00 = ExpHash(rhash(x0, y0, seed1), rhash(x0, y0, seed2), max, omega)
    m01 = ExpHash(rhash(x0, y0 + 1, seed1), rhash(x0, y0 + 1, seed2), max, omega)
    m10 = ExpHash(rhash(x0 + 1, y0, seed1), rhash(x0 + 1, y0, seed2), max, omega)
    m11 = ExpHash(rhash(x0 + 1, y0 + 1, seed1), rhash(x0 + 1, y0 + 1, seed2), max, omega)

    uax = igt[:uax]
    vax = igt[:vax]
    ubx = igt[:ubx]
    vbx = igt[:vbx]
    uay = igt[:uay]
    vay = igt[:vay]
    uby = igt[:uby]
    vby = igt[:vby]

    FillUp!(uax, m00 * cos(b00), n)
    FillDn!(vax, m01 * cos(b01), n)

    FillUp!(ubx, m10 * cos(b10), n)
    FillDn!(vbx, m11 * cos(b11), n)

    FillUp!(uay, m00 * sin(b00), n)
    FillUp!(vay, m01 * sin(b01), n)

    FillDn!(uby, m10 * sin(b10), n)
    FillDn!(vby, m11 * sin(b11), n)

    # FillUp!(uax, 1.0 * cos(b00), n)
    # FillDn!(vax, 1.0 * cos(b01), n)
    #
    # FillUp!(ubx, 1.0 * cos(b10), n)
    # FillDn!(vbx, 1.0 * cos(b11), n)
    #
    # FillUp!(uay, 1.0 * sin(b00), n)
    # FillUp!(vay, 1.0 * sin(b01), n)
    #
    # FillDn!(uby, 1.0 * sin(b10), n)
    # FillDn!(vby, 1.0 * sin(b11), n)

    rv = @dict uax vax ubx vbx uay vay uby vby

    return rv
end

function init_edge_tables_amn(x0::Int64, y0::Int64, n::Int64, omega::Float64, seed1, seed2, igt)
    b00 = rhash(x0, y0, seed1)
    b01 = rhash(x0, y0 + 1, seed1)
    b10 = rhash(x0 + 1, y0, seed1)
    b11 = rhash(x0 + 1, y0 + 1, seed1)

    max = typemax(UInt64)

    # m00 = ExpHash(rhash(x0, y0, seed1), rhash(x0, y0, seed2), max, omega)
    # m01 = ExpHash(rhash(x0, y0 + 1, seed1), rhash(x0, y0 + 1, seed2), max, omega)
    # m10 = ExpHash(rhash(x0 + 1, y0, seed1), rhash(x0 + 1, y0, seed2), max, omega)
    # m11 = ExpHash(rhash(x0 + 1, y0 + 1, seed1), rhash(x0 + 1, y0 + 1, seed2), max, omega)

    uax = igt[:uax]
    vax = igt[:vax]
    ubx = igt[:ubx]
    vbx = igt[:vbx]
    uay = igt[:uay]
    vay = igt[:vay]
    uby = igt[:uby]
    vby = igt[:vby]

    # FillUp!(uax, m00 * cos(b00), n)
    # FillDn!(vax, m01 * cos(b01), n)
    #
    # FillUp!(ubx, m10 * cos(b10), n)
    # FillDn!(vbx, m11 * cos(b11), n)
    #
    # FillUp!(uay, m00 * sin(b00), n)
    # FillUp!(vay, m01 * sin(b01), n)
    #
    # FillDn!(uby, m10 * sin(b10), n)
    # FillDn!(vby, m11 * sin(b11), n)

    FillUp!(uax, 1.0 * cos(b00), n)
    FillDn!(vax, 1.0 * cos(b01), n)

    FillUp!(ubx, 1.0 * cos(b10), n)
    FillDn!(vbx, 1.0 * cos(b11), n)

    FillUp!(uay, 1.0 * sin(b00), n)
    FillUp!(vay, 1.0 * sin(b01), n)

    FillDn!(uby, 1.0 * sin(b10), n)
    FillDn!(vby, 1.0 * sin(b11), n)

    rv = @dict uax vax ubx vbx uay vay uby vby

    return rv
end

function get_noise(i::Int64, j::Int64, igt, spline)
    uax = igt[:uax]
    vax = igt[:vax]
    ubx = igt[:ubx]
    vbx = igt[:vbx]
    uay = igt[:uay]
    vay = igt[:vay]
    uby = igt[:uby]
    vby = igt[:vby]

    u = uax[j] + uay[i]
    v = vax[j] + vay[i]

    a = lerp(spline[j], u, v)

    w = ubx[j] + uby[i]
    z = vbx[j] + vby[i]

    b = lerp(spline[j], w, z)

    return lerp(spline[i], a, b)
end

function get_noise(n::Int64, i0::Int64, j0::Int64, cell, igt, spline)
    for i = 1:n
        for j = 1:n
            cell[i0+i, j0+j] = get_noise(i, j, igt, spline)
        end
    end
end

function add_noise(n::Int64, i0::Int64, j0::Int64, scale::Float64, cell, igt, spline)
    for i = 1:n
        for j = 1:n
            cell[i0+i, j0+j] += scale * get_noise(i, j, igt, spline)
        end
    end
end

function Amortised_Exp_Noise_2D_cell(x::Int64, y::Int64, n::Int64, seed_1::Int64, seed_2::Int64, omega::Float64, cell, scale::Float64, m0::Int64, m1::Int64)
    #
    #
    # s_curve(x) = x^2 * (3 - 2 * x)
    # lerp(t, a, b) = a + t * (b - a)
    #
    # B = 0x100
    # BM = 0xff
    # N = 0x1000

    uax = Array{Float64,1}(undef, n)
    vax = Array{Float64,1}(undef, n)

    ubx = Array{Float64,1}(undef, n)
    vbx = Array{Float64,1}(undef, n)

    uay = Array{Float64,1}(undef, n)
    vay = Array{Float64,1}(undef, n)

    uby = Array{Float64,1}(undef, n)
    vby = Array{Float64,1}(undef, n)

    igt = @dict uax vax ubx vbx uay vay uby vby

    r = 1

    for i = 1:(m0-1)
        n = fld(n, 2)
        r += r
    end

    sp_tab = init_spline_table(n)

    for i0 = 0:(r-1)
        for j0 = 0:(r-1)
            igt = init_edge_tables_exp(x + i0, y + j0, n, omega, seed_1, seed_2, igt)
            get_noise(n, i0 * n, j0 * n, cell, igt, sp_tab)
        end
    end

    k = m0
    while (k < m1 && n >= 2)
        n = fld(n, 2)
        r += r
        x += x
        y += y
        spline_table_inplace!(sp_tab, n)
        for i0 = 0:(r-1)
            for j0 = 0:(r-1)
                igt = init_edge_tables_exp(x + i0, y + j0, n, omega, seed_1, seed_2, igt)
                add_noise(n, i0 * n, j0 * n, scale, cell, igt, sp_tab)
            end
        end
        scale *= scale
        k = k + 1
    end
end

function Amortised_Noise_2D_cell(x::Int64, y::Int64, n::Int64, seed_1::Int64, seed_2::Int64, omega::Float64, cell, scale::Float64, m0::Int64, m1::Int64)
    #
    #
    # s_curve(x) = x^2 * (3 - 2 * x)
    # lerp(t, a, b) = a + t * (b - a)
    #
    # B = 0x100
    # BM = 0xff
    # N = 0x1000

    uax = Array{Float64,1}(undef, n)
    vax = Array{Float64,1}(undef, n)

    ubx = Array{Float64,1}(undef, n)
    vbx = Array{Float64,1}(undef, n)

    uay = Array{Float64,1}(undef, n)
    vay = Array{Float64,1}(undef, n)

    uby = Array{Float64,1}(undef, n)
    vby = Array{Float64,1}(undef, n)

    igt = @dict uax vax ubx vbx uay vay uby vby

    r = 1

    for i = 1:(m0-1)
        n = fld(n, 2)
        r += r
    end

    sp_tab = init_spline_table(n)

    for i0 = 0:(r-1)
        for j0 = 0:(r-1)
            igt = init_edge_tables_amn(x + i0, y + j0, n, omega, seed_1, seed_2, igt)
            get_noise(n, i0 * n, j0 * n, cell, igt, sp_tab)
        end
    end

    k = m0
    while (k < m1 && n >= 2)
        n = fld(n, 2)
        r += r
        x += x
        y += y
        spline_table_inplace!(sp_tab, n)
        for i0 = 0:(r-1)
            for j0 = 0:(r-1)
                igt = init_edge_tables_amn(x + i0, y + j0, n, omega, seed_1, seed_2, igt)
                add_noise(n, i0 * n, j0 * n, scale, cell, igt, sp_tab)
            end
        end
        scale *= scale
        k = k + 1
    end
end

function AMN_EXP_2D_grid(coords, seed_1, seed_2, omega, scale, m0, m1, cellsize)
    cell = zeros(cellsize, cellsize)

    Amortised_Exp_Noise_2D_cell(coords[1], coords[2], cellsize, seed_1, seed_2, omega, cell, scale, m0, m1)
    return cell
end

function AMN_2D_grid(coords, seed_1, seed_2, omega, scale, m0, m1, cellsize)
    cell = zeros(cellsize, cellsize)

    Amortised_Noise_2D_cell(coords[1], coords[2], cellsize, seed_1, seed_2, omega, cell, scale, m0, m1)
    return cell
end

# grd = AMN_EXP_2D_grid((12, 11), 12, 456768, 0.1, 0.65, 1, 5, 128)
# # grd = grd[1:1023, 1:1023]
# plot(1:128, 1:128, grd, st = :heatmap, c = :grays)
