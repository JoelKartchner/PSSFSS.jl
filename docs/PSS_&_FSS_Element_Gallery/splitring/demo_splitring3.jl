# ---
# title: Square Split Ring with Narrow Gap
# cover: "assets/demo_splitring3.png"
# description: "Square narrow-gap split rings in a square lattice created by the splitring function"
# ---

using PSSFSS, Plots

units = mm
s1 = [4, 0]
s2 = [0, 4]
gapwidth = [0.6]
gapcenter = [0]
sides = 4
a = [1.626]
b = [2.475]
sheet = splitring(; units, sides, orient = 45, a, b, ntri=400, s1, s2, gapwidth, gapcenter)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:red, size=(400,400), rep=(4,4)) #src
savefig("assets/demo_splitring3.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))
