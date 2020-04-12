using SharedArrays
using Distributed
using Random
using Plots
using JLD
using StatsBase


d = load("$(ARGS[1]).jld")
NE0D = d["NE0D"]

labels = ["50/50" "Mostly Negative" "Mostly Positive"]
colors = ["blue" "red" "green"]
pyplot()

for i = 1:3
    prop = collect(sort(proportionmap(NE0D[i,:])))
    xs = [pair[1] for pair in prop]
    ys = [pair[2] * 100 for pair in prop]
    mask = xs .> 1
    xs = xs[mask]
    ys = ys[mask]
    plot!(xs, ys, label=labels[i], yaxis = :log, seriescolor = colors[i])
    plot!(xs, ys, seriestype=:scatter, yaxis=:log, seriescolor = colors[i], label="")
end

plot!(title = "Julia Validation with ns=$(size(NE0D)[2])", dpi=300)
plot!(xlabel = "Richness", ylabel = "Frequency (%)")
xticks!(2:maximum(NE0D))
yticks!([0.001, 0.01, 0.1, 1, 10, 100])
ylims!((0.001, 100.))

png("$(ARGS[1])")
