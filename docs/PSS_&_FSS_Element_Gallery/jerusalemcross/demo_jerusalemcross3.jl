# ---
# title: Partially Loaded Jerusalem Cross
# cover: "assets/demo_jerusalemcross3.png"
# description: "Partially loaded jerusalem cross in a square lattice, created by the jerusalemcross function"
# ---

using PSSFSS, Plots
sheet = jerusalemcross(P=1, L1=0.9, L2=0.12, w=0.04, A = 0.4, B = 0.08, units=cm, ntri=800) 
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:green, size=(400,400), rep=(4,4)) #src
savefig("assets/demo_jerusalemcross3.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))
