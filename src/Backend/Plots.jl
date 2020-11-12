module Plots

using DataFrames, Plots

default(dpi = 300, titlefont = (16, "times"), guidefont = (12, "times"), tickfont = (12, "times"), legendfont = (6, "times"))

function runchart(raw::DataFrame, data::Symbol, batch::Symbol; 
                center::Function = mean, 
                spec::Real = 100,
                USL::Real = 110,
                LSL::Real = 90,
                sigma::Real = 3,
                title::String = "Run chart",
                xlabel::String = "Batch",
                ylabel::String = "Assay (%)")
    
    gdf = groupby(raw, batch)
    point = length(unique(raw[!, batch]))     
    run = size(raw, 1) / point
    plot([1, point], repeat([USL], 2),color = :green, linestyle = :dash, label = "USL/LSL", legend = :topright)
    plot!([1, point] ,repeat([LSL], 2),color = :green, linestyle = :dash, label = false)
    plot!([1, point], repeat([median(combine(gdf, data => center)[!,2])], 2), color = :black, label = "Median", linestyle = :dash )
    if sigma > 0
        plot!([1, point], repeat([spec + sigma * std(raw[!, data])], 2), color = :chartreuse, linestyle = :dash, label = "±$(sigma)σ")
        plot!([1, point], repeat([spec - sigma * std(raw[!, data])], 2), color = :chartreuse, linestyle = :dash, label = false)
    end
    scatter!(raw[!, batch], raw[!, data], label = false, color = :deepskyblue, markeralpha = .3)
    plot!(1:point, combine(gdf, data => center)[!,2], color = :black, width = 3, label = "Batch $center")
    up = max(USL, raw[!, data]...)
    down = min(LSL, raw[!, data]...)
    Δ = up - down
    μ = (up + down) / 2
    ylims!(μ - 1.5 * Δ, μ + 2 * Δ)
    xlabel!(xlabel)
    ylabel!(ylabel)
    title!(title)
end

    
function runchart(raw::DataFrame, cpkresult::Vextor{Real}, data::Symbol, batch::Symbol; 
                center::Function = mean, 
                spec::Real = 100,
                USL::Real = 110,
                LSL::Real = 90,
                sigma::Real = 3,
                title::String = "Run chart",
                xlabel::String = "Batch",
                ylabel::String = "Assay (%)")

    gdf = groupby(raw, batch)
    point = length(unique(raw[!, batch]))     
    run = size(raw, 1) / point
    plot([1, point], repeat([USL], 2),color = :green, linestyle = :dash, label = "USL/LSL", legend = :topright)
    plot!([1, point] ,repeat([LSL], 2),color = :green, linestyle = :dash, label = false)
    plot!([1, point], repeat([median(combine(gdf, data => center)[!,2])], 2), color = :black, label = "Median", linestyle = :dash )
    if sigma > 0
    plot!([1, point], repeat([spec + sigma * std(raw[!, data])], 2), color = :chartreuse, linestyle = :dash, label = "±$(sigma)σ")
    plot!([1, point], repeat([spec - sigma * std(raw[!, data])], 2), color = :chartreuse, linestyle = :dash, label = false)
    end
    good = findall(x->x > 1.33, cpkresult)
    medium = findall(x->(x < 1.33) & (x > 1), cpkresult)
    bad = findall(x->x < 1, cpkresult)
    isempty(good) || begin
        sub = filter(batch => x->x in good, raw)
        scatter!(sub[!, batch], sub[!, data], label = "Cpk > 1.33", color = :deepskyblue, markeralpha = .3)
    end
    isempty(medium) || begin
        sub = filter(batch => x->x in medium, raw)
        scatter!(sub[!, batch], sub[!, data], label = "1 < Cpk < 1.33", color = :red, markeralpha = .3)
    end
    isempty(bad) || begin
        sub = filter(batch => x->x in bad, raw)
        scatter!(sub[!, batch], sub[!, data], label = "Cpk < 1", color = :purple, markeralpha = .3)
    end
    plot!(1:point, combine(gdf, data => center)[!,2], color = :black, width = 3, label = "Batch $center")
    up = max(USL, raw[!, data]...)
    down = min(LSL, raw[!, data]...)
    Δ = up - down
    μ = (up + down) / 2
    ylims!(μ - 1.5 * Δ, μ + 2 * Δ)
    xlabel!(xlabel)
    ylabel!(ylabel)
    title!(title)
end


# xbarchart
# ichart: n(subgroup) = 1

# rchart: n(subgroup) = 2-8
# schart: n(subgroup) > 8
# mrchart: n(subgroup) = 1

# capability histogram
# capability plot


end