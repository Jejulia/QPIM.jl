module QPIM

export plot_init

include("Backend\\Stats.jl")

# Plot module
# Avoid import Plots as default
const path = @__DIR__()
plot_init() = include("$path\\Backend\\Plots.jl")

__init__() = precompile(plot_init,())

end
