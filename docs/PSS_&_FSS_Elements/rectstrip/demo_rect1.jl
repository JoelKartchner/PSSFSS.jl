# ---
# title: Square Patch
# cover: "assets/demo_rect1.png"
# description: "Square patch created by the rectstrip function"
# ---

# The `rectstrip` element employs a structured mesh, so 
# it can be analyzed very efficiently.

using Plots, PSSFSS
patch = rectstrip(Nx=10, Ny=10, Px=1, Py=1, Lx=0.5, Ly=0.5, units=cm)
p1 = plot(patch, unitcell=true, linecolor=:red)
plot(patch, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:green, size=(400,300), rep=(4,3)) #src
savefig("assets/demo_rect1.png") #src
p2 = plot(patch, rep=(3,3), linecolor=:blue)
plot(p1, p2, layout=(1,2))
