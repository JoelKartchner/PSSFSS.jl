# ---
# title: Hexagonal Loop with Patch
# cover: "assets/demo_polyring4.png"
# description: "Hexagonal loop with patch in a triangular lattice, created by the polyring function"
# ---

using PSSFSS, Plots
sheet = polyring(s1=[1, 0], s2=[0.5, √3/2], a=[0.0, 0.45], b=[0.15, 0.55], sides=6, units=cm, orient=30, ntri=600)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:orange, size=(400,400), rep=(4,4)) #src
savefig("assets/demo_polyring4.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))
