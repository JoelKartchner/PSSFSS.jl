# ---
# title: Horizontal Grid
# cover: "assets/demo_rect2.png"
# description: "Horizontal grid created by the rectstrip function"
# ---

# The [`rectstrip`](@ref) element employs a structured mesh, so 
# it can be analyzed very efficiently.

using Plots, PSSFSS
strip = rectstrip(Nx=4, Ny=8, Px=0.2, Py=1, Lx=0.2, Ly=0.15, units=cm)
p1 = plot(strip, unitcell=true, linecolor=:red)
plot(strip, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:green, size=(400,300), rep=(15,3)) #src
savefig("assets/demo_rect2.png") #src
p2 = plot(strip, rep=(5,3), linecolor=:blue)
plot(p1, p2, layout=(1,2))
