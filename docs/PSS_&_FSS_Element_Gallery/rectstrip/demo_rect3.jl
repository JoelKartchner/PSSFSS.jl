# ---
# title: Vertical Grid
# cover: "assets/demo_rect3.png"
# description: "Vertical grid created by the rectstrip function"
# ---

# The [`rectstrip`](@ref) element employs a structured mesh, so 
# it can be analyzed very efficiently.

using Plots, PSSFSS
strip = rectstrip(Nx=8, Ny=4, Px=1, Py=0.2, Lx=0.15, Ly=0.2, units=cm)
p1 = plot(strip, unitcell=true, linecolor=:red)
plot(strip, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:green, size=(400,300), rep=(4,15)) #src
savefig("assets/demo_rect3.png") #src
p2 = plot(strip, rep=(3,5), linecolor=:blue)
plot(p1, p2, layout=(1,2), size=(800,400))
