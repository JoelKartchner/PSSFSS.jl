# ---
# title: Circular Split Rings with Wide Gaps
# cover: "assets/demo_splitring2.png"
# description: "Circular wide-gap split rings in a square lattice created by the splitring function"
# ---

using PSSFSS, Plots

a = [4.7, 6.8]
b = [5.9, 8]
s1 = [19.5, 0]
s2 = [0, 19.5]
gapangle = [26, 85.6]
gapcenter = [90, 0]
sheet = splitring(; units=mm, sides=42, ntri=800, a, b, s1, s2, gapangle, gapcenter)
p1 = plot(sheet, linecolor=:red, unitcell=true)
p2 = plot(sheet, linecolor=:blue, rep=(4,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:red, size=(400,400), rep=(4,4)) #src
savefig("assets/demo_splitring2.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))
