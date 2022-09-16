# ---
# title: Diagonal Strip in Square Lattice
# cover: "assets/demo_diagstrip1.png"
# description: "Diagonally oriented strip in a square lattice, created by the diagstrip function"
# ---

# The [`diagstrip`](@ref) element requires a square unit cell and restricts the values of `orient`
# to either `45` or `-45`.  Using this element (rather than a rotated `rectstrip`) is useful when
# there are other sheets in the FSS/PSS structure that all share the same square lattice.  
# In this case, the interactions between
# sheets can be rigorously accounted for using higher-order generalized scattering parameters.

using PSSFSS, Plots
sheet = diagstrip(P=5.2, w=0.21, units=mm, Nl=60, Nw=4, orient=45)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:red, size=(400,400), rep=(4,4)) #src
savefig("assets/demo_diagstrip1.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))
