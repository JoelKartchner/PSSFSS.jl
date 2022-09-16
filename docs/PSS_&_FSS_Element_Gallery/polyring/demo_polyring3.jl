# ---
# title: Circular Loop
# cover: "assets/demo_polyring3.png"
# description: "Circular loop in a triangular lattice, created by the polyring function"
# ---

# We approximate a circular ring as a 40-sided polygon, placed in a close-packed, triangular lattice.

using PSSFSS, Plots
sheet = polyring(s1=[1, 0], s2=[0.5, √3/2], a=[0.25], b=[0.4], sides=40, units=cm, ntri=600)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:orange, size=(400,400), rep=(4,4)) #src
savefig("assets/demo_polyring3.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))
