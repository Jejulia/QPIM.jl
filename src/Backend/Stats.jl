module Stats

# Backend module for statistics functionality

using Statistics, Distributions, DataFrames

export 

    runaboutmedian, eaboutmedian, runupdown, eupdown,

    cluster, mixture, trend, oscillation,

    cpk, ppk

# ------------------------------------------------------------------------------------------------
# Run chart analysis
function runaboutmedian(points::Vector{Real})
    n = 1
    center = median(points)
    for (i,v) in enumerate(points[2:end])
       (v - center) * (points[i] - center) < 0 && (n += 1)
    end
    n
 end

function eaboutmedian(points::Vector{Real})
    N = length(points)
    m = length(findall(x->x > median(points), points))
    n = N - m
    1 + 2 * m * n / N
end

function cluster(points::Vector{Real})
    N = length(points)
    R = runaboutmedian(points)
    E = eaboutmedian(points)
    Z = (R - E) / sqrt((E - 1) * (E - 2) / (N - 1))
    cdf(Normal(0, 1), Z)
end

mixture(data::Vector{Real}) = 1 - cluster(data)

function runupdown(points::Vector{Real})
    n = 1
    df = diff(points)
    prod = df[1] > 0
    for i in df[2:end]
        i  != 0 && (prod * i < 0) && (prod *= -1 ; n += 1)
    end
    n
end

eupdown(points::Vector{Real}) = 2 / 3 * (2 * length(points) - 1)

function trend(points::Vector{Real})
    N = length(points)
    v = runupdown(points)
    Z = (v - (2 * N - 1) / 3) / sqrt((16 * N  - 29) / 90)
    cdf(Normal(0, 1), Z)
end

oscillation(points::Vector{Real}) = 1 - trend(points)

# ---------------------------------------------------------------------------------------------------------
# Process capability/performance index
# cp/pp/z.bench/z.usl/z.lsl/ppm

function cpk(x::Vector{Real}; s::Real = 3, USL::Real = 110, LSL::Real = 90, correct::Bool = true)
    d = length(x) - 1
    μ = mean(x)
    correct ? (sd = std(x) / (sqrt(2 / d) * gamma(d / 2 + 0.5) / gamma(d / 2))) : (sd = std(x))
    min((USL - μ) / s / sd, (μ - LSL) / s / sd)
 end
 
function cpk(x::Matrix{Real}; s::Real = 3 USL::Real = 110, LSL::Real = 90, correct::Bool = true)
    sd = map(1:size(x, 2)) do i
       std(x[:,i])
    end
    d = size(x, 1) * size(x, 2) - size(x, 2)
    μ = mean(x)
    sd = sum(sd .* (size(x, 1) - 1)) / d
    correct && (sd /= (sqrt(2 / d) * gamma(d / 2 + 0.5) / gamma(d / 2)))
    min((USL - μ) / s / sd, (μ - LSL) / s  / sd)
end

cpk(x::DataFrame; s::Real = 3 USL::Real = 110, LSL::Real = 90, correct::Bool = true) = 
    cpk(Matrix{Real}(x); s = s, USL = USL, LSL = LSL, correct = correct)

 
function ppk(x::Matrix{Real}; s::Real = 3, USL::Real = 110, LSL::Real = 90, correct::Bool = false)
    d = size(x, 1) * size(x, 2) - 1
    μ = mean(x)
    correct ? (sd = std(x) / (sqrt(2 / d) * gamma(d / 2 + 0.5) / gamma(d / 2))) : (sd = std(x))
    min((USL - μ) / s / sd, (μ - LSL) / s / sd)
end

ppk(x::DataFrame; s::Real, USL::Real = 110, LSL::Real = 90, correct::Bool = false) = 
    ppk(Matrix{Real}(x), s = s, USL = USL, LSL = LSL, correct = correct)

end