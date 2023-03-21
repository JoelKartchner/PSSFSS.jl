# ---
# title: Sinuous (3-Arm)
# cover: "assets/demo_sinuous1.png"
# description: "Sinuous 3-arm element in a hexagonal lattice, created by the sinuous function"
# ---

using PSSFSS, Plots
P = 0.75
s1 = P * [1, 0]; s2 = P * [0.5, √3/2]; units = cm
b = [0.1, 0.18,  0.267,  0.35]; orient = 90
sheet = sinuous(; arms=3, b, orient, w=0.03, rc=0.05, g=0.035, sides=45, ntri=1500, units, s1, s2)

p1 = plot(sheet, linecolor=:red, size=(400,400), unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(3,3))
plot(sheet, axis=false, xlabel="", ylabel="", grid=false, linecolor=:blue, size=(400,400), rep=(3,3)) #src
savefig("assets/demo_sinuous1.png") #src
plot(p1, p2, layout = (1,2), size=(800,400), margin=10Plots.pt)
