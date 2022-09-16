# ---
# title: Loaded Cross
# cover: "assets/demo_loadedcross2.png"
# description: "Loaded Cross in a square lattice, created by the loadedcross function"
# ---

using PSSFSS, Plots
sheet = loadedcross(s1=[1,0], s2=[0,1], L1=0.9, L2=0.25, w=0.06, units=cm, ntri=400)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:blue, size=(400,400), rep=(4,4)) #src
savefig("assets/demo_loadedcross2.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))
