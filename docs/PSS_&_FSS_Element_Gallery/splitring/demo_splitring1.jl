# ---
# title: Circular Split Rings with Narrow Double Gaps
# cover: "assets/demo_splitring1.png"
# description: "Circular split rings in a square lattice created by the splitring function"
# ---

using PSSFSS, Plots

units = mm
a = [1.9, 5.7]
b = a + [1.4, 1.3]
gapwidth = [(1.5, 1.5), (1.5, 1.5)]
gapcenter = [(45,225), (45,225)]
s1 = [15, 0]
s2 = [0, 15]
sheet = splitring(; units=mm, sides=45, ntri=700, a, b, s1, s2, gapwidth, gapcenter)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:red, size=(400,400), rep=(4,4)) #src
savefig("assets/demo_splitring1.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))
