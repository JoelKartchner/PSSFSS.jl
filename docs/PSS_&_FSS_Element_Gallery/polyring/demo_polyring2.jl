# ---
# title: Octagonal Loop in Octagonal Hole
# cover: "assets/demo_polyring2.png"
# description: "Octagonal loop within an octagonal hole, created by the polyring function"
# ---

# Here we have an octagonal hole and an octagonal, annular loop.  
# Such a geometry could be more efficiently modeled
# by making use of a sheet of class `M`, and triangulating the complementary
# (smaller in area) region to that shown below.  However, 
# if one needs to model this type of inductive element while 
# including the surface resistance of the sheet, 
# this can only be done using a sheet of class `J` (the default class).

using PSSFSS, Plots
sheet = polyring(s1=[1,0], s2=[0,1], a=[0.2, 0.5], b=[0.35, -25], sides=8, orient=22.5, units=cm, ntri=400)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:orange, size=(400,400), rep=(4,3)) #src
savefig("assets/demo_polyring2.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))
