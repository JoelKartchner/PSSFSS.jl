# ---
# title: Meanderline Polarizer Sheet
# cover: "assets/demo_meander1.png"
# description: "Meanderline trace created by the meander function"
# ---

# Like the [`rectstrip`](@ref) element, the [`meander`](@ref) element employs a structured mesh, so 
# it can be analyzed very efficiently.

using PSSFSS, Plots
sheet = meander(a=0.3535, b=0.707, h=0.28, w1=0.018, w2=0.018, ntri=600, units=inch,rot=45)
p1 = plot(sheet, linecolor=:green, unitcell=true)
p2 = plot(sheet, linecolor=:green, rep=(4,3))
plot(sheet, axis=false, xlabel="", ylabel="", xtick=[], ytick=[], linecolor=:red, size=(400,400), rep=(4,3)) #src
savefig("assets/demo_meander1.png") #src
plot(p1, p2, layout = (1,2), size=(800,400))
