```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/examples.jl"
```

# [PSSFSS](https://github.com/simonp0420/PSSFSS) Examples

```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/symmetric_strip.jl"
```

## Symmetric Strip Grating
This example consists of a symmetric strip grating, i.e. a grating where the strip width
is half the unit cell period ``P``:

![diagram](./assets/symmetric_strip_diagram.png)

Only three of the infinite number of strips in the grating are shown, and they extend infinitely to the left and right.
The grating lies in the ``z=0`` plane with free space on both sides. The shaded areas represent metalization.
The dashed lines show two possible choices for the unit cell location: "J" for a formulation in terms of electric
surface currents, and "M" for magnetic surface currents.

For normal incidence there is a closed-form solution due to Weinstein,
but for a more recent reference one can find the solution in Problem 10.6 of R. E. Collin,
*Field Theory of Guided Waves, Second Ed.*,
IEEE Press, 1991.  Here is the code for computing the exact solution:

````@example symmetric_strip
"""
    grating(kP, nterms=30) -> (Γ, T)

Compute the normal incidence reflecton and transmission coefficients of a symmetric grid of
zero-thickness conducting strips.  The product of the period of the strips and the incident
electric field wavenumber is `kP` (dimensionless).  The incident electric field is
perpendicular to the direction along the axis of the strips.  The series have been
accelerated by applying a Kummer's transformation, using the first two terms in the Maclaurin
series for the inverse sin function.  `kP` must be in the half-open interval [0,1). The
default number of summed terms `nterms` yields better than 10 digits of accuracy over the
interval [0.01,0.99].
"""
function grating(kP; nterms=30)
    sum1 = 1.3862943611198906 # \sum_{n=1}^{\infty} 1/(n-1/2) - 1/n = log(4)
    sum3 = 7.2123414189575710 # \sum_{n=1}^{\infty} (n-1/2)^{-3} - n^{-3} = 6 * \zeta(3)
    x = kP/(4π)
    θ = x*sum1 + x^3/6 * sum3
    for n = 1:nterms
        xonmhalf = x/(n - 0.5)
        xon = x/n
        term = asin(xonmhalf) - (xonmhalf + (xonmhalf)^3/6) -
              (asin(xon) - (xon + xon^3/6))
        θ += term
    end
    Γ = sin(θ) * cis(-π/2 - θ)
    T = 1 + Γ
    return (Γ, T)
end
````

 Note that using the extension of
[Babinet's Principle for electromagnetic fields](http://kirkmcd.princeton.edu/examples/babinet.pdf)
this also provides the solution (upon appropriate interchange and sign change of the coefficients) for
the case where the incident wave polarization is parallel to the direction of the strips.

Here is the PSSFSS code to analyze this structure using electric currents as the unknowns.  We will
scale the geometry so that the frequency in GHz is numerically equal to the period of the strips
measured in wavelengths.

````@example symmetric_strip
using Plots, PSSFSS
c = 11.802852677165355 # light speed [inch*GHz]
period = c  # so the period/wavelength = freq in GHz
Py = period
Ly = period/2
Px = Lx = Ly/10 # Infinite in x direction so this can be anything
Ny = 60
Nx = round(Int, Ny*Lx/Ly)
sheet = rectstrip(;Px, Py, Lx, Ly, Nx, Ny, units=inch)
flist = 0.02:0.02:0.98
steering = (θ=0, ϕ=0)
strata = [Layer()
          sheet
          Layer()]
results_j = analyze(strata, flist, steering, showprogress=false,
                    resultfile=devnull, logfile=devnull);
p1 = plot(sheet)
p2 = plot(sheet, unitcell=true)
ptitle = plot(title = "Symmetric Strip Triangulation",
             grid = false, showaxis = false, xtick=[], ytick=[],
             bottom_margin = -50Plots.px)
plot(ptitle, p1, p2, layout = @layout([A{0.09h}; [B C]]))
savefig("symstrip1.png"); nothing  # hide
````

![](symstrip1.png)

Note that setting `Lx = Px` causes the strip to fully occupy the x-extent
of the unit cell.  PSSFSS automatically ensures that the triangle edges at these unit
cell boundaries define basis functions that satisfy the Floquet (phase shift) boundary
conditions, so that currents are free to flow across these unit cell boundaries.

We can also analyze the same structure using magnetic currents in the areas free of
metalization as the unknowns:

````@example symmetric_strip
sheet = rectstrip(;class='M', Px, Py, Lx, Ly, Nx, Ny, units=inch)
strata = [Layer()
          sheet
          Layer()]
results_m = analyze(strata, flist, steering, showprogress=false,
                    resultfile=devnull, logfile=devnull);
nothing #hide
````

Each 50-frequency run of `analyze` takes about 14 seconds
for this geometry of 720 triangles on my machine.
More detailed timing information is available in the log file
(which is omitted for generating this documentation).

We will compare the PSSFSS results to the analytic solution:

````@example symmetric_strip
# Generate exact results:
rt = grating.(2π*flist)
rperp_exact = first.(rt)
tperp_exact = last.(rt)
rpar_exact = -tperp_exact
tpar_exact = -rperp_exact;
nothing #hide
````

Obtain PSSFSS results for electric and magnetic currents:

````@example symmetric_strip
outrequest = @outputs s11(v,v) s21(v,v) s11(h,h) s21(h,h)
rperp_j, tperp_j, rpar_j, tpar_j =
      collect.(eachcol(extract_result(results_j, outrequest)))
rperp_m, tperp_m, rpar_m, tpar_m =
      collect.(eachcol(extract_result(results_m, outrequest)));
nothing #hide
````

Generate the comparison plots:

````@example symmetric_strip
angdeg(z) = rad2deg(angle(z)) # Convenience function

p1 = plot(title = "Perpendicular Reflection Magnitude",
          xlabel = "Period (wavelengths)",
          ylabel = "Coefficient Magnitude",
          legend=:topleft)
plot!(p1, flist, abs.(rperp_exact), ls=:dash, label="Exact")
plot!(p1, flist, abs.(rperp_j), label="PSSFSS J")
plot!(p1, flist, abs.(rperp_m), label="PSSFSS M")
savefig("symstrip2.png"); nothing  # hide
````

![](symstrip2.png)

````@example symmetric_strip
p2 = plot(title = "Perpendicular Reflection Phase",
          xlabel = "Period (wavelengths)",
          ylabel = "Phase (deg)")
plot!(p2, flist, angdeg.(rperp_exact), ls=:dash, label="Exact")
plot!(p2, flist, angdeg.(rperp_j), label="PSSFSS J")
plot!(p2, flist, angdeg.(rperp_m), label="PSSFSS M")
savefig("symstrip3.png"); nothing  # hide
````

![](symstrip3.png)

````@example symmetric_strip
p1 = plot(title = "Parallel Reflection Magnitude",
          xlabel = "Period (wavelengths)",
          ylabel = "Coefficient Magnitude")
plot!(p1, flist, abs.(rpar_exact), ls=:dash, label="Exact")
plot!(p1, flist, abs.(rpar_j), label="PSSFSS J")
plot!(p1, flist, abs.(rpar_m), label="PSSFSS M")
savefig("symstrip4.png"); nothing  # hide
````

![](symstrip4.png)

````@example symmetric_strip
p2 = plot(title = "Parallel Reflection Phase",
          xlabel = "Period (wavelengths)",
          ylabel = "Phase (deg)")
plot!(p2, flist, angdeg.(rpar_exact), ls=:dash, label="Exact")
plot!(p2, flist, angdeg.(rpar_j), label="PSSFSS J")
plot!(p2, flist, angdeg.(rpar_m), label="PSSFSS M")
savefig("symstrip5.png"); nothing  # hide
````

![](symstrip5.png)

Now look at the transmission coefficients:

````@example symmetric_strip
p1 = plot(title = "Perpendicular Transmission Magnitude",
          xlabel = "Period (wavelengths)",
          ylabel = "Coefficient Magnitude")
plot!(p1, flist, abs.(tperp_exact), ls=:dash, label="Exact")
plot!(p1, flist, abs.(tperp_j), label="PSSFSS J")
plot!(p1, flist, abs.(tperp_m), label="PSSFSS M")
savefig("symstrip6.png"); nothing  # hide
````

![](symstrip6.png)

````@example symmetric_strip
p2 = plot(title = "Perpendicular Transmission Phase",
          xlabel = "Period (wavelengths)",
          ylabel = "Phase (deg)")
plot!(p2, flist, angdeg.(tperp_exact), ls=:dash, label="Exact")
plot!(p2, flist, angdeg.(tperp_j), label="PSSFSS J")
plot!(p2, flist, angdeg.(tperp_m), label="PSSFSS M")
savefig("symstrip7.png"); nothing  # hide
````

![](symstrip7.png)

````@example symmetric_strip
p1 = plot(title = "Parallel Transmission Magnitude",
          xlabel = "Period (wavelengths)",
          ylabel = "Coefficient Magnitude", legend=:topleft)
plot!(p1, flist, abs.(tpar_exact), ls=:dash, label="Exact")
plot!(p1, flist, abs.(tpar_j), label="PSSFSS J")
plot!(p1, flist, abs.(tpar_m), label="PSSFSS M")
savefig("symstrip8.png"); nothing  # hide
````

![](symstrip8.png)

````@example symmetric_strip
p2 = plot(title = "Parallel Transmission Phase",
          xlabel = "Period (wavelengths)",
          ylabel = "Phase (deg)")
plot!(p2, flist, angdeg.(tpar_exact), ls=:dash, label="Exact")
plot!(p2, flist, angdeg.(tpar_j), label="PSSFSS J")
plot!(p2, flist, angdeg.(tpar_m), label="PSSFSS M")
savefig("symstrip9.png"); nothing  # hide
````

![](symstrip9.png)

### Conclusion
Although very good agreement is obtained, as expected the best agreement between
all three results occurs for the lowest frequencies, where the triangles are
smallest in terms of wavelength.  This emphasizes the fact that it is necessary for the
user to check that enough triangles have been requested for good convergence
over the frequency band of interest.  This example is an extremely demanding case
in terms of bandwidth, as the ratio of maximum to minimum frequency here
is ``0.98/0.02 = 49:1``

```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/resistive_square_patch.jl"
```

## Resistive Square Patch
This example will demonstrate the ability of PSSFSS to accurately model finite
conductivity of FSS metalization.  It consists of a square finitely conducting
patch in a square lattice.  It is taken from a paper by Alon S. Barlevy and
Yahya Rahmat-Samii,
"Fundamental Constraints on the Electrical Characteristics of Frequency Selective
Surfaces", **Electromagnetics**, vol. 17, 1997, pp. 41-68. This particular example
is from Section 3.2, Figures 7 and 8.  We will compare PSSFSS results to those digitized
from the cited figures.

We start by defining a function that creates a patch of the desired sheet resistance:

````@example resistive_square_patch
using Plots, PSSFSS
patch(R) = rectstrip(Nx=10, Ny=10, Px=1, Py=1, Lx=0.5, Ly=0.5, units=cm, Rsheet=R)
plot(patch(0), unitcell=true)
savefig("resistive1.png"); nothing  # hide
````

![](resistive1.png)

The patches measure 0.5 cm on a side and lie in a square lattice of period 1 cm.
Now we perform the analysis, looping over the desired values of sheet resistance.

````@example resistive_square_patch
steering = (ϕ=0, θ=0)
flist = 1:0.5:60
Rs = [0, 10, 30, 100]
calculated = zeros(length(flist), length(Rs)) # preallocate storage
outputs = @outputs s11mag(v,v)
for (i,R) in pairs(Rs)
    strata = [Layer(), patch(R), Layer()]
    results = analyze(strata, flist, steering, showprogress=false,
                      logfile=devnull, resultfile=devnull)
    calculated[:,i] = extract_result(results, outputs)
end
````

Looping over the four sheet resistance values, each evaluated at 119 frequencies
required approximately 20 seconds on my machine.

We plot the results, including those digitized from the paper for comparison:

````@example resistive_square_patch
using DelimitedFiles
markers = (:diamond, :utriangle, :square, :xcross)
colors = (:blue, :red, :green, :black)
p = plot(xlim=(-0.01,60.01), xtick = 0:10:60, ylim=(-0.01,1.01), ytick=0:0.1:1,
         xlabel="Frequency (GHz)", ylabel="Reflection Coefficient Magnitude",
         title = "Resistive Square Patch",
         legend=:topright)
for (i,R) in pairs(Rs)
    scatter!(p, flist, calculated[:,i], label="PSSFSS $R Ω", ms=2, shape=markers[i], color=colors[i])
    data = readdlm("../src/assets/barlevy_patch_digitized_$(R)ohm.csv", ',')
    plot!(p, data[:,1], data[:,2], label="Barlevy $R Ω", ls=:dash, color=colors[i])
end
p
savefig("resistive2.png"); nothing  # hide
````

![](resistive2.png)

### Conclusion
PSSFSS results are indistinguishable from those reported in the cited paper.

```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/cross_on_dielectric_substrate.jl"
```

## Cross on Dielectric Substrate
This example is also taken from the paper by Alon S. Barlevy and
Yahya Rahmat-Samii, "Fundamental Constraints on the Electrical Characteristics
of Frequency Selective Surfaces", **Electromagnetics**, vol. 17, 1997, pp. 41-68.
This particular example is from Section 3.2, Figures 7 and 8.  It also appeared at
higher resolution in Barlevy's PhD dissertation from which the comparison curves
were digitized.

We use the `loadedcross` element where we choose `w > L2/2`, so that the Cross
is "unloaded", i.e. the center section is filled in with metalization:

````@example cross_on_dielectric_substrate
using Plots, PSSFSS, DelimitedFiles
sheet = loadedcross(w=1.0, L1=0.6875, L2=0.0625, s1=[1.0,0.0],
                    s2=[0.0,1.0], ntri=600, units=cm)
plot(sheet, unitcell=true)
savefig("cross1.png"); nothing  # hide
````

![](cross1.png)

A few things to note. First, the mesh is **unstructured**.  So there are no redundant
triangle face-pairs that PSSFSS can exploit to reduce execution time.  Second, the
number of triangle faces generated is only approximately equal to the requested value
of 600.  This can be verified by entering the Julia variable `sheet` at the
[REPL](https://docs.julialang.org/en/v1/manual/getting-started/#man-getting-started)
(i.e. the Julia prompt):

````@example cross_on_dielectric_substrate
sheet
````

Alternatively, the `facecount` function will return the number of triangle faces on the sheet:

````@example cross_on_dielectric_substrate
facecount(sheet)
````

The cross FSS is etched on a dielectric sheet of thickness 3 mm.  The dielectric
constant is varied over the values 1, 2, and 4 to observe the effect on the resonant
frequency.  Following the reference, the list of analysis frequencies is varied slightly
depending on the value of dielectric constant:

````@example cross_on_dielectric_substrate
resultsstack = Any[]
steering = (ϕ=0, θ=0)
for eps in [1, 2, 4]
    strata = [  Layer()
                sheet
                Layer(ϵᵣ=eps, width=3mm)
                Layer()
             ]
    if eps == 1
        flist = 1:0.2:30
    elseif eps == 2
        flist = 1:0.2:26
    else
        flist = 1:0.2:20
    end
    results = analyze(strata, flist, steering, showprogress=false,
                      resultfile=devnull, logfile=devnull)
    push!(resultsstack, results)
end
````

The above loop requires about 85 seconds of execution time on my machine.
Compare PSSFSS results to those digitized from the dissertation figure:

````@example cross_on_dielectric_substrate
col=[:red,:blue,:green]
p = plot(xlim=(0.,30), xtick = 0:5:30, ylim=(0,1), ytick=0:0.1:1,
         xlabel="Frequency (GHz)", ylabel="Reflection Coefficient Magnitude",
         legend=:topleft, lw=2)
for (i,eps) in enumerate([1,2,4])
    data = extract_result(resultsstack[i],  @outputs FGHz s11mag(v,v))
    plot!(p, data[:,1], data[:,2], label="PSSFSS ϵᵣ = $eps", lc=col[i])
    data = readdlm("../src/assets/barlevy_diss_eps$(eps).csv", ',')
    plot!(p, data[:,1], data[:,2], label="Barlevy ϵᵣ = $eps", lc=col[i], ls=:dot)
end
p
savefig("cross2.png"); nothing  # hide
````

![](cross2.png)

### Conclusion
PSSFSS results agree very well with those of the cited reference, especially when
accounting for the fact that the reference results are obtained by digitizing a
scanned figure.

```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/square_loop_absorber.jl"
```

## Square Loop Absorber
This example is from Figure 7 of Costa and Monorchio: "A frequency selective
radome with wideband absorbing properties", *IEEE Trans. AP-S*,
Vol. 60, no. 6, June 2012, pp. 2740--2747.  It shows how one can use the `polyring`
function to model square loop elements.  Three different designs are examined
that employ different loop thicknesses and different values of sheet resistance.
We compare the reflection coefficient magnitudes computed by PSSFSS with those digitized
from the cited figure when the sheet is suspended
5 mm above a ground plane, hence we will also make use of the `pecsheet` function.

````@example square_loop_absorber
using Plots, PSSFSS, DelimitedFiles
D = 11 # Period of square lattice (mm)
r_outer = √2/2 * D/8 * [5,6,7] # radii of square outer vertices
thickness = D/16 * [1,2,3]
r_inner = r_outer - √2 * thickness  # radii of square inner vertices
Rs = [15,40,70] # Sheet resistance (Ω/□)
labels = ["Thin", "Medium", "Thick"]
colors = [:green, :blue, :red]
p = plot(title="Costa Absorber", xlim=(0,25),ylim=(-35,0),xtick=0:5:25,ytick=-35:5:0,
         xlabel="Frequency (GHz)", ylabel="Reflection Magnitude (dB)", legend=:bottomleft)
ps = []
for (i,(ri, ro, label, color, R)) in enumerate(zip(r_inner, r_outer, labels, colors, Rs))
    sheet = polyring(sides=4, s1=[D, 0], s2=[0, D], ntri=700, orient=45,
                     a=[ri], b=[ro], Rsheet=R, units=mm)
    push!(ps, plot(sheet, unitcell=true, title=label, lc=color))
    strata = [Layer()
              sheet
              Layer(width=5mm)
              pecsheet() # Perfectly conducting ground plane
              Layer()]
    results = analyze(strata, 1:0.2:25, (ϕ=0, θ=0), showprogress=false,
                      resultfile=devnull, logfile=devnull)
    data = extract_result(results, @outputs FGHz s11dB(h,h))
    plot!(p, data[:,1], data[:,2], label="PSSFSS "*label, lc=color)
    dat = readdlm("../src/assets/costa_2014_" * lowercase(label) * "_reflection.csv", ',')
    plot!(p, dat[:,1], dat[:,2], label="Costa "*label, ls=:dash, lc=color)
end
plot(ps..., layout=(1,3))
savefig("sqloop1.png"); nothing  # hide
````

![](sqloop1.png)

This run takes about 83 seconds on my machine.

````@example square_loop_absorber
p
savefig(p,"sqloop2.png"); nothing  # hide
````

![](sqloop2.png)

It is useful to take a look at the log file created by PSSFSS for the last run above:
```
Starting PSSFSS 1.0.0 analysis on 2022-09-14 at 14:31:14.261
Julia Version 1.8.1
Commit afb6c60d69a (2022-09-06 15:09 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 8 × Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.1 (ORCJIT, skylake)
  Threads: 8 on 8 virtual cores
  BLAS: LBTConfig([ILP64] libopenblas64_.so)



Dielectric layer information...

 Layer  Width  units  epsr   tandel   mur  mtandel modes  beta1x  beta1y  beta2x  beta2y
 ----- ------------- ------- ------ ------- ------ ----- ------- ------- ------- -------
     1    0.0000  mm    1.00 0.0000    1.00 0.0000     2   571.2    -0.0    -0.0   571.2
 ==================  Sheet   1  ========================   571.2    -0.0    -0.0   571.2
     2    5.0000  mm    1.00 0.0000    1.00 0.0000    42   571.2    -0.0    -0.0   571.2
 ==================  Sheet   2  ========================     0.0     0.0     0.0     0.0
     3    0.0000  mm    1.00 0.0000    1.00 0.0000     2   571.2    -0.0    -0.0   571.2



PSS/FSS sheet information...

Sheet  Loc         Style      Rot  J/M Faces Edges Nodes Unknowns  NUFP
-----  ---  ---------------- ----- --- ----- ----- ----- -------- ------
   1     1          polyring   0.0  J    753  1201   448    1058  567009
   2     2              NULL   0.0  E      0     0     0       0       0

...
```

Note from the dielectric layer report that there are 42 modes defined in the region between the
ground plane and the FSS sheet.  This is the number of modes selected by the code to include
in the generalized scattering matrix formulation to properly account for electromagnetic coupling
between the two surfaces. If the 5 mm spacing were increased to, say, 7 mm then fewer modes
would be needed.  Also note in the FSS sheet information that `NUFP` (the number of unique face pairs)
is exactly equal to the number of faces squared (``567009 = 753^2``), a consequence of the unstructured
triangulation used for a `polyring`.

### Conclusion
PSSFSS results agree very well with those of the paper, except for the medium
width loop, where the agreement is not quite as good.  The reason for this is
not known.

```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/splitringexample.jl"
```

## Split-Ring Resonator
This example is taken from Figure 3 of
Fabian‐Gongora, et al, "Independently Tunable Closely Spaced Triband
Frequency Selective Surface Unit Cell Using the Third Resonant Mode of Split Ring Slots",
IEEE Access, 8/3/2021, Digital Object Identifier 10.1109/ACCESS.2021.3100325.

It consists of three concentric split rings with gaps sequentially rotated by 180° situated
on a thin dielectric slab.

Here is a script that analyzes this geometry:

```Julia
using PSSFSS, Plots, DelimitedFiles

r123 = [3.7, 4.25, 4.8]
d = 10.8
w = 0.3
g = 0.3
a = r123 .- w/2
b = a .+ w
gapwidth = [g, g, g]
gapcenter = [-90, 90, -90]
s1 = [d, 0]
s2 = [0, d]

sheet = splitring(;class='M', units=mm, sides=42, ntri=1500, a, b, s1, s2, gapwidth, gapcenter)
display(plot(sheet, unitcell=true))
strata = [
    Layer()
    sheet
    Layer(ϵᵣ=10.2, tanδ=0.0023, width=0.13mm)
    Layer()
    ]

freqs = union(1:0.1:2, 2.02:0.02:3, 3.1:0.05:14)
steering = (θ = 0, ϕ = 0)
results = analyze(strata, freqs, steering)
s11dbvv = extract_result(results, @outputs s11db(v,v))
p = plot(xlabel="Frequency (GHz)", ylabel="S₁₁ Amplitude (dB)",
    xtick=1:14, xlim=(1, 14), ytick = -40:5:0, ylim=(-30,0),
    legend=:bottom)
plot!(p, freqs, s11dbvv, label="PSSFSS")
dat = readdlm("../src/assets/fabian2021_fig3_digitized.csv", ',')
plot!(p, dat[:,1], dat[:,2], label="Fabian (CST)")
```

![](./assets/fabian2021_element.png)

![](./assets/fabian2021_comparison.png)

Generally good agreement is seen between the PSSFSS predicted reflection amplitude and that
digitized from the paper. (The latter was obtained from a CST frequency domain analysis,
according to the paper's authors.) However, there is a small discrepancy in the predicted resonant
frequencies that increases
with frequency, likely because both results are less well converged at higher frequencies.

```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/band_pass_filter.jl"
```

## Loaded Cross Band Pass Filter
This example is originally from Fig. 7.9 of B. Munk, *Frequency Selective Surfaces,
Theory and Design,* John Wiley and Sons, 2000.  The same case was analyzed in L. Li,
D. H. Werner et al, "A Model-Based Parameter Estimation Technique for
Wide-Band Interpolation of Periodic Moment Method Impedance Matrices With Application to
Genetic Algorithm Optimization of Frequency Selective Surfaces", *IEEE Trans. AP-S*,
vol. 54, no. 3, March 2006, pp. 908-924, Fig. 6.  Unfortunately, in neither reference
are the dimensions of the loaded cross stated, except for the square unit cell
period of 8.61 mm.  I estimated the dimensions from the sketch in Fig. 6 of the second
reference.  To provide a reliable comparison, I enlisted my colleague
[Mike Maybell](https://www.linkedin.com/in/mike-maybell-308b77ba),
principal of Planet Earth Communications, who generously offered to
analyze the filter using
[CST Microwave Studio](https://www.3ds.com/products-services/simulia/products/cst-studio-suite/),
a rigorous commercial finite volume electromagnetic solver.

Two identical loaded cross slot-type elements are separated by a 6 mm layer of dielectric
constant 1.9.  Outboard of each sheet is a 1.1 cm layer of dielectric constant 1.3.
The closely spaced sheets are a good test of the generalized scattering formulation
implemented in PSSFSS.  The sheet geometry is shown below.  Remember that the entire
sheet is metalized *except* for the region of the triangulation.

````@example band_pass_filter
using Plots, PSSFSS
sheet = loadedcross(class='M', w=0.023, L1=0.8, L2=0.14,
            s1=[0.861,0.0], s2=[0.0,0.861], ntri=600, units=cm)
plot(sheet, unitcell=true)
savefig("bpf1.png"); nothing  # hide
````

![](bpf1.png)

````@example band_pass_filter
steering = (ϕ=0, θ=0)
strata = [  Layer()
            Layer(ϵᵣ=1.3, width=1.1cm)
            sheet
            Layer(ϵᵣ=1.9, width=0.6cm)
            sheet
            Layer(ϵᵣ=1.3, width=1.1cm)
            Layer()  ]
flist = 1:0.1:20
results = analyze(strata, flist, steering, resultfile=devnull,
                  logfile=devnull, showprogress=false)
data = extract_result(results, @outputs FGHz s21db(v,v) s11db(v,v))
using DelimitedFiles
dat = readdlm("../src/assets/MaybellLoadedCrossResults.csv", ',', skipstart=1)
p = plot(xlabel="Frequency (GHz)", ylabel="Reflection Coefficient (dB)",
         legend=:left, title="Loaded Cross Band-Pass Filter", xtick=0:2:20, ytick=-30:5:0,
         xlim=(-0.1,20.1), ylim=(-35,0.1))
plot!(p, data[:,1], data[:,3], label="PSSFSS", color=:red)
plot!(p, dat[:,1], dat[:,2], label="CST", color=:blue)
savefig("bpf2.png"); nothing  # hide
````

![](bpf2.png)

````@example band_pass_filter
p2 = plot(xlabel="Frequency (GHz)", ylabel="Transmission Coefficient (dB)",
          legend=:bottom, title="Loaded Cross Band-Pass Filter", xtick=0:2:20, ytick=-80:10:0,
         xlim=(-0.1,20.1), ylim=(-80,0.1))
plot!(p2, data[:,1], data[:,2], label="PSSFSS", color=:red)
plot!(p2, dat[:,1], dat[:,4], label="CST", color=:blue)
savefig("bpf3.png"); nothing  # hide
````

![](bpf3.png)

This analysis takes about 57 seconds for 191 frequencies on my machine.  Note that
rather than including two separate invocations of the `loadedcross` function when
defining the strata, I referenced the same sheet object in the two different locations.
This allows PSSFSS to recognize that the triangulations are identical, and to exploit
this fact in making the analysis more efficient.  In fact, if both sheets had been embedded
in similar dielectric claddings (in the same order), then the GSM (generalized scattering matrix)
computed for the sheet in its first location could be reused without additional computation for its
second location.  In this case, though, only the spatial integrals are re-used.  For an oblique
incidence case, computing the spatial integrals is often the most expensive part of the analysis,
so the savings from reusing the same sheet definition can be substantial.

### Conclusion
Very good agreement is obtained versus CST over a large dynamic range.

```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/cpss_optimization.jl"
```

## Meanderline-Based CPSS
A "CPSS" is a circular polarization selective structure, i.e., a structure that reacts differently
to the two senses of circular polarization.
We'll first look at analyzing a design presented in the literature, and then proceed to optimize another
design using PSSFSS as the analysis engine inside the optimization objective function.
### Sjöberg and Ericsson Design
This example comes from the paper D. Sjöberg and A. Ericsson, "A multi layer meander line circular
polarization selective structure (MLML-CPSS)," The 8th European Conference on Antennas and Propagation
(EuCAP 2014), 2014, pp. 464-468, doi: 10.1109/EuCAP.2014.6901792.

The authors describe an ingenious structure consisting of 5 progressively rotated meanderline sheets, which
acts as a circular polarization selective surface: it passes LHCP (almost) without attenuation or
reflection, and reflects RHCP (without changing its sense!) almost without attenuation or transmission.

Here is the script that analyzes their design:

````@example cpss_optimization
using PSSFSS
# Define convenience functions for sheets:
outer(rot) = meander(a=3.97, b=3.97, w1=0.13, w2=0.13, h=2.53+0.13, units=mm, ntri=600, rot=rot)
inner(rot) = meander(a=3.97*√2, b=3.97/√2, w1=0.1, w2=0.1, h=0.14+0.1, units=mm, ntri=600, rot=rot)
center(rot) = meander(a=3.97, b=3.97, w1=0.34, w2=0.34, h=2.51+0.34, units=mm, ntri=600, rot=rot)
# Note our definition of `h` differs from that of the reference by the width of the strip.
t1 = 4mm # Outer layers thickness
t2 = 2.45mm # Inner layers thickness
substrate = Layer(width=0.1mm, epsr=2.6)
foam(w) = Layer(width=w, epsr=1.05) # Foam layer convenience function
rot0 = 0 # rotation of first sheet
strata = [
    Layer()
    outer(rot0)
    substrate
    foam(t1)
    inner(rot0 - 45)
    substrate
    foam(t2)
    center(rot0 - 2*45)
    substrate
    foam(t2)
    inner(rot0 - 3*45)
    substrate
    foam(t1)
    outer(rot0 - 4*45)
    substrate
    Layer() ]
steering = (θ=0, ϕ=0)
flist = 10:0.1:20
#
results = analyze(strata, flist, steering, showprogress=false,
                  resultfile=devnull, logfile=devnull);
nothing #hide
````

Here are plots of the five meanderline sheets:

````@example cpss_optimization
using Plots
plot(outer(rot0), unitcell=true, title="Sheet1")
````

````@example cpss_optimization
plot(inner(rot0-45), unitcell=true, title="Sheet2")
````

````@example cpss_optimization
plot(center(rot0-2*45), unitcell=true, title="Sheet3 (Center)")
````

````@example cpss_optimization
plot(inner(rot0-3*45), unitcell=true, title="Sheet4")
````

````@example cpss_optimization
plot(outer(rot0-4*45), unitcell=true, title="Sheet5")
````

Notice that not only are the meanders rotated, but so too are the unit cell rectangles.
This is because we used the generic `rot` keyword argument that rotates the entire unit
cell and its contents. `rot` can be used for any FSS or PSS element type.  As a consequence of
the different rotations applied to each unit cell, interactions between sheets due to higher-order
modes cannot be accounted for; only the dominant ``m=n=0`` TE and TM modes are used in cascading
the individual sheet scattering matrices.  This approximation is adequate for sheets that are
sufficiently separated.  We can see from the log file (saved from a previous run where it was not
disabled) that only 2 modes are used to model the interactions between sheets:

```
Starting PSSFSS 1.0.0 analysis on 2022-09-12 at 16:35:20.914
Julia Version 1.8.1
Commit afb6c60d69 (2022-09-06 15:09 UTC)
Platform Info:
  OS: Windows (x86_64-w64-mingw32)
  CPU: 8 × Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.1 (ORCJIT, skylake)
  Threads: 8 on 8 virtual cores
  BLAS: LBTConfig([ILP64] libopenblas64_.dll)



******************* Warning ***********************
   Unequal unit cells in sheets 1 and 2
   Setting #modes in dividing layer 3 to 2
******************* Warning ***********************


******************* Warning ***********************
   Unequal unit cells in sheets 2 and 3
   Setting #modes in dividing layer 5 to 2
******************* Warning ***********************


******************* Warning ***********************
   Unequal unit cells in sheets 3 and 4
   Setting #modes in dividing layer 7 to 2
******************* Warning ***********************


******************* Warning ***********************
   Unequal unit cells in sheets 4 and 5
   Setting #modes in dividing layer 9 to 2
******************* Warning ***********************


Dielectric layer information...

 Layer  Width  units  epsr   tandel   mur  mtandel modes  beta1x  beta1y  beta2x  beta2y
 ----- ------------- ------- ------ ------- ------ ----- ------- ------- ------- -------
     1    0.0000  mm    1.00 0.0000    1.00 0.0000     2  1582.7    -0.0    -0.0  1582.7
 ==================  Sheet   1  ========================  1582.7    -0.0    -0.0  1582.7
     2    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
     3    4.0000  mm    1.05 0.0000    1.00 0.0000     2  1582.7    -0.0    -0.0  1582.7
 ==================  Sheet   2  ========================   791.3  -791.3  1582.7  1582.7
     4    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
     5    2.4500  mm    1.05 0.0000    1.00 0.0000     2   791.3  -791.3  1582.7  1582.7
 ==================  Sheet   3  ========================     0.0 -1582.7  1582.7     0.0
     6    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
     7    2.4500  mm    1.05 0.0000    1.00 0.0000     2     0.0 -1582.7  1582.7     0.0
 ==================  Sheet   4  ========================  -791.3  -791.3  1582.7 -1582.7
     8    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
     9    4.0000  mm    1.05 0.0000    1.00 0.0000     2  -791.3  -791.3  1582.7 -1582.7
 ==================  Sheet   5  ======================== -1582.7     0.0     0.0 -1582.7
    10    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
    11    0.0000  mm    1.00 0.0000    1.00 0.0000     2 -1582.7     0.0     0.0 -1582.7

...
```

Note that PSSFSS prints warnings to the log file where it is forced to set the number of layer
modes to 2 because of unequal unit cells.  Also, in the dielectric layer list it can be seen
that these layers are assigned 2 modes each.  The thin layers adjacent to sheets are assigned 0
modes because these sheets are incorporated into so-called "GSM blocks" or "Gblocks" wherein
the presence of the thin layer is accounted for using the stratified medium Green's functions.

Here is the script that compares PSSFSS predicted performance with very
high accuracy predictions from CST and COMSOL that were digitized from figures in the paper.

````@example cpss_optimization
using Plots, DelimitedFiles
RL11rr = -extract_result(results, @outputs s11db(r,r))
AR11r = extract_result(results, @outputs ar11db(r))
IL21L = -extract_result(results, @outputs s21db(L,L))
AR21L = extract_result(results, @outputs ar21db(L))

RLgoal, ILgoal, ARgoal  = ([0.5, 0.5], [0.5, 0.5], [0.75, 0.75])
foptlimits = [12, 18]
default(lw=2, xtick=10:20, xlabel="Frequency (GHz)", ylabel="Amplitude (dB)", gridalpha=0.3)

p = plot(flist,RL11rr,title="RHCP → RHCP Return Loss", label="PSSFSS",
         ylim=(0,2), ytick=0:0.25:2)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_rl.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_rl.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, RLgoal, color=:black, lw=4, label="Goal")
savefig("cpssa1.png"); nothing  # hide
````

![](cpssa1.png)

````@example cpss_optimization
p = plot(flist,AR11r,title="RHCP → RHCP Reflected Axial Ratio", label="PSSFSS",
         xlim=(10,20), ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_ar_reflected.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_ar_reflected.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, ARgoal, color=:black, lw=4, label="Goal")
savefig("cpssa2.png"); nothing  # hide
````

![](cpssa2.png)

````@example cpss_optimization
p = plot(flist,IL21L,title="LHCP → LHCP Insertion Loss", label="PSSFSS",
         xlim=(10,20), ylim=(0,2), ytick=0:0.25:2)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_il.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_il.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, ILgoal, color=:black, lw=4, label="Goal")
savefig("cpssa3.png"); nothing  # hide
````

![](cpssa3.png)

````@example cpss_optimization
p = plot(flist,AR21L,title="LHCP → LHCP Transmitted Axial Ratio", label="PSSFSS",
         xlim=(10,20), ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_ar_transmitted.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_ar_transmitted.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, ARgoal, color=:black, lw=4, label="Goal")
savefig("cpssa4.png"); nothing  # hide
````

![](cpssa4.png)

The PSSFSS results generally track well with the high-accuracy solutions, but are less accurate
especially at the high end of the band, possibly because in PSSFSS metallization thickness is
neglected and cascading is performed for this structure using only the two principal Floquet modes.
As previosly discussed, this is necessary because the rotated meanderlines are achieved by rotating
the entire unit cell, and the unit cell for sheets 2 and 4 are not square.  Since the periodicity of
the sheets in the structure varies from sheet to sheet, higher order Floquet modes common to neighboring
sheets cannot be defined, so we are forced to use only the dominant (0,0) modes which are independent of
the periodicity.  This limitation is removed in a later example.
Meanwhile, it is of interest to note that their high-accuracy runs
required 10 hours for CST and 19 hours for COMSOL on large engineering workstations.  The PSSFSS run
took about 50 seconds on my desktop machine.

### Design Based on PSSFSS Optimization with CMAES
Here we use PSSFSS in conjunction with the CMAES optimizer from the
[CMAEvolutionStrategy](https://github.com/jbrea/CMAEvolutionStrategy.jl) package.  I've used CMAES
in the past with good success on some tough optimization problems.  Here is the code that defines
the objective function:

```julia
using PSSFSS
using Dates: now

let bestf = typemax(Float64)
    global objective
    """
        result = objective(x)

    """
    function objective(x)
        ao, bo, ai, bi, ac, bc, wo, ho, wi, hi, wc, hc, t1, t2 = x
        (bo > ho > 2.1*wo && bi > hi > 2.1*wi && bc > hc > 2.1*wc) || (return 5000.0)

        outer(rot) = meander(a=ao, b=bo, w1=wo, w2=wo, h=ho, units=mm, ntri=400, rot=rot)
        inner(rot) = meander(a=ai, b=bi, w1=wi, w2=wi, h=hi, units=mm, ntri=400, rot=rot)
        center(rot) = meander(a=ac, b=bc, w1=wc, w2=wc, h=hc, units=mm, ntri=400, rot=rot)

        substrate = Layer(width=0.1mm, epsr=2.6)
        foam(w) = Layer(width=w, epsr=1.05)
        rot0 = 0

        strata = [
                Layer()
                outer(rot0)
                substrate
                foam(t1*1mm)
                inner(rot0 - 45)
                substrate
                foam(t2*1mm)
                center(rot0 - 2*45)
                substrate
                foam(t2*1mm)
                inner(rot0 - 3*45)
                substrate
                foam(t1*1mm)
                outer(rot0 - 4*45)
                substrate
                Layer() ]
        steering = (θ=0, ϕ=0)
        flist = 11.5:0.25:18.5
        resultfile = logfile = devnull
        showprogress = false
        results = analyze(strata, flist, steering; showprogress, resultfile, logfile)
        s11rr, s21ll, ar11db, ar21db = eachcol(extract_result(results,
                       @outputs s11db(R,R) s21db(L,L) ar11db(R) ar21db(L)))
        RL = -s11rr
        IL = -s21ll
        RLgoal, ILgoal, ARgoal  = (0.4, 0.5, 0.6)
        obj = max(maximum(RL) - RLgoal,
                  maximum(IL) - ILgoal,
                  maximum(ar11db) - ARgoal,
                  maximum(ar21db) - ARgoal)
        if obj < bestf
            bestf = obj
            open("optimization_best.log", "a") do fid
                xround = map(t  -> round(t, digits=4), x)
                println(fid, round(obj,digits=5), " at x = ", xround, "  #", now())
            end
        end
        return obj
    end
end
```

We optimize at 29 frequencies between 11.5 and 18.5 GHz.  In the previously presented
Sjöberg and Ericsson design a smaller frequency range of 12 to 18 GHz was used for optimization.
We have also adopted more ambitious goals for return loss, insertion loss, and axial ratio
of 0.4 dB, 0.5 dB, and 0.6 dB, respectively. These should be feasible because for
our optimization setup, we have relaxed the restriction that all unit cells must be
identical squares.  With more "knobs to turn" (i.e., a larger search space), we expect to be able
to do somewhat better in terms of bandwidth and worst-case performance.  For our setup there are
fourteen optimization variables in all.
As one can see from the code above, each successive sheet in the structure is rotated an additional
45 degrees relative to its predecessor.   The objective is defined as the maximum departure
from the goal of RHCP reflected return loss, LHCP insertion loss, or reflected or transmitted axial ratio that
occurs at any of the analysis frequencies (i.e. we are setting up for "minimax" optimization). Also,
the `let` block allows the objective function to maintain persistent state in the
variable `bestf` which is initialized to the largest 64-bit floating point value. Each time a set
of inputs results in a lower objective function value, `bestf` is updated with this value and
the inputs and objective function value are
written to the file "optimization_best.log", along with a date/time stamp.  This allows the user
to monitor the optimization and to terminate the it prematurely, if desired, without losing the
best result achieved so far. Each objective function evaluation takes about 9 seconds on my machine.

Here is the code for running the optimization:

```julia
using CMAEvolutionStrategy
#  x = [a1,  b1,   a2, b2,  a3,  b3,  wo,  ho,  wi,   hi,  wc,   hc,  t1,  t2]
xmin = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.1, 0.1, 0.1,  0.1, 0.1,  0.1, 2.0, 2.0]
xmax = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 0.35,4.0, 0.35, 4.0, 0.35, 4.0, 6.0, 6.0]
x0 = 0.5 * (xmin + xmax)
isfile("optimization_best.log") && rm("optimization_best.log")
popsize = 2*(4 + floor(Int, 3*log(length(x0))))
opt = minimize(objective, x0, 1.0;
           lower = xmin,
           upper = xmax,
           maxfevals = 9000,
           popsize = popsize,
           xtol = 1e-4,
           ftol = 1e-6)
```

Note that I set the population size to twice the normal default value.  Based
on previous experience, using 2 to 3 times the default population size helps the
optimizer to do better on tough objective functions like the present one.
I let the optimizer run overnight for about 17 hours, during which time it reduced the objective function
value from 11.48 dB to -0.12 dB. Since it appeared to have essentially converged, I terminated the run.

Here is a look at the beginning and final portions of the file "optimization_best.log":
```
11.47502 at x = [4.3371, 3.6377, 4.6978, 3.0065, 4.2178, 3.7151, 0.3456, 1.4562, 0.3345, 2.1578, 0.3498, 2.0437, 3.3593, 3.5733]  #2022-09-13T14:35:39.320
8.83597 at x = [4.3381, 4.8048, 3.4625, 4.7309, 3.051, 3.5356, 0.1072, 3.585, 0.1039, 2.2732, 0.3206, 1.1203, 3.8363, 2.2236]  #2022-09-13T14:35:57.550
6.61653 at x = [3.4206, 3.1036, 3.0215, 3.2064, 3.0201, 3.4882, 0.1034, 1.5586, 0.103, 2.3935, 0.2026, 2.1608, 4.2708, 2.4485]  #2022-09-13T14:36:06.768
6.17284 at x = [4.6927, 3.9841, 4.6078, 3.8153, 3.2699, 3.0023, 0.3254, 3.0281, 0.3167, 1.0311, 0.3334, 1.3243, 4.9482, 2.5076]  #2022-09-13T14:37:01.227
4.14975 at x = [4.4617, 4.105, 3.9944, 4.8857, 4.6484, 3.6307, 0.3457, 2.0488, 0.3354, 1.9338, 0.3275, 1.5484, 3.2037, 2.683]  #2022-09-13T14:37:17.647
2.2661 at x = [4.7303, 4.1483, 3.4428, 3.3463, 3.607, 3.2735, 0.1572, 2.2432, 0.35, 2.2474, 0.2468, 1.376, 2.6822, 2.3709]  #2022-09-13T14:39:58.103
1.57851 at x = [4.95, 4.1857, 4.0174, 3.3486, 4.9004, 3.678, 0.12, 1.6691, 0.2451, 1.0522, 0.2964, 1.7673, 3.6622, 2.6612]  #2022-09-13T14:41:33.641
...
-0.11824 at x = [3.0642, 4.8748, 4.1393, 3.0066, 3.8585, 3.0047, 0.344, 3.0002, 0.3469, 1.4042, 0.2502, 2.3751, 3.766, 2.3372]  #2022-09-14T06:37:45.783
-0.1184 at x = [3.0416, 4.8731, 4.148, 3.0012, 3.8553, 3.0103, 0.3449, 3.0104, 0.3482, 1.4088, 0.2533, 2.3858, 3.7593, 2.3347]  #2022-09-14T06:38:57.587
-0.11896 at x = [3.0454, 4.881, 4.1304, 3.0001, 3.8389, 3.01, 0.345, 3.0042, 0.348, 1.413, 0.2524, 2.3719, 3.7682, 2.3398]  #2022-09-14T06:50:21.349
-0.11961 at x = [3.069, 4.8773, 4.1652, 3.002, 3.8446, 3.0055, 0.3439, 3.0101, 0.3468, 1.4051, 0.2518, 2.3774, 3.7652, 2.3401]  #2022-09-14T06:59:09.133
-0.11961 at x = [3.0661, 4.8788, 4.155, 3.0014, 3.8314, 3.0069, 0.3441, 2.9933, 0.3476, 1.4168, 0.2516, 2.3689, 3.7816, 2.3322]  #2022-09-14T07:06:13.332
-0.1197 at x = [3.0657, 4.8757, 4.1549, 3.0018, 3.8331, 3.0061, 0.3436, 3.0024, 0.3468, 1.4055, 0.2517, 2.3716, 3.7758, 2.3382]  #2022-09-14T07:15:07.056
```

The final sheet geometries and performance of this design are shown below:

![](./assets/cpss_cmaesopt_sheets.png)

![](./assets/cpss_cmaesopt_rl_refl.png)

![](./assets/cpss_cmaesopt_ar_refl.png)

![](./assets/cpss_cmaesopt_il_trans.png)

![](./assets/cpss_cmaesopt_ar_trans.png)

As hoped for, the performance meets the more stringent design goals over a broader bandwidth than the
Sjöberg and Ericsson design, presumably because of the greater design flexibility allowed here.

```@meta
EditURL = "https://github.com/simonp0420/PSSFSS.jl/tree/main/docs/literate/cpss2.jl"
```

## Meanderline/Strip-Based CPSS
This example comes from the same authors as the previous example.  The paper is
A. Ericsson and D. Sjöberg, "Design and Analysis of a Multilayer Meander Line
Circular Polarization Selective Structure", IEEE Trans. Antennas Propagat.,
Vol. 65, No. 8, Aug 2017, pp. 4089-4101.
The design is similar to that of the previous example except that here, the two ``\pm 45^\circ``
rotated meanderlines are replaced with rectangular strips.
This allows us to employ the `diagstrip` element and the `orient` keyword for the
`meander` elements to maintain the same, square unit cell for all sheets. By doing this
we allow PSSFSS to rigorously account for the inter-sheet coupling using multiple
high-order modes in the generalized scattering matrix (GSM) formulation.

We begin by computing the skin depth and sheet resistance for the
copper traces.  The conductivity and thickness are as stated in the paper:

````@example cpss2
# Compute skin depth and sheet resistance:
using PSSFSS.Constants: μ₀ # free-space permeability [H/m]
f = (10:0.1:20) * 1e9 # frequencies in Hz
σ = 58e6 # conductivity of metalization [S/m]
t = 18e-6 # metalization thickness [m]
Δ = sqrt.(2 ./ (2π*f*σ*μ₀)) # skin depth [m]
@show extrema(t./Δ)
````

````@example cpss2
Rs = 1 ./ (σ * Δ)
@show extrema(Rs)
````

We see that the metal is many skin depths thick (effectively infinitely thick) so that we can use
the thick metal surface sheet resistance formula.  Since the latter varies with frequency, we approximate
it over the band 10-20 GHz by a value near its mean: 0.032 Ω/□.

Here is the script that analyzes the design from the referenced paper:

````@example cpss2
using PSSFSS
P = 5.2 # side length of unit square
d1 = 2.61 # Inner layer thickness
d2 = 3.81 # Outer layer thickness
h0 = 2.44 # Inner meanderline dimension (using paper's definition of h)
h2 = 2.83 # Outer meanderline dimension (using paper's definition of h)
w0x = 0.46 # Inner meanderline line thickness of traces running along x
w0y = 0.58 # Inner meanderline line thickness of traces running along y
w1 = 0.21 # Rectangular strips width
w2x = 0.25   # Outer meanderline line thickness of traces running along x
w2y = 0.17 # Outer meanderline line thickness of traces running along y

outer(orient) = meander(a=P, b=P, w1=w2y, w2=w2x, h=h2+w2x, units=mm, ntri=600,
                        Rsheet=0.032, orient=orient)
inner = meander(a=P, b=P, w1=w0y, w2=w0x, h=h0+w0x, units=mm, ntri=600, Rsheet=0.032)
strip(orient) = diagstrip(P=P, w=w1, units=mm, Nl=60, Nw=4, orient=orient, Rsheet=0.032)

substrate = Layer(width=0.127mm, epsr=2.17, tandel=0.0009)
foam(w) = Layer(width=w, epsr=1.043, tandel=0.0017)
sheets = [outer(-90), strip(-45), inner, strip(45), outer(90)]
strata = [
    Layer()
    substrate
    sheets[1]
    foam(d2 * 1mm)
    substrate
    sheets[2]
    foam(d1 * 1mm)
    sheets[3]
    substrate
    foam(d1 * 1mm)
    substrate
    sheets[4]
    foam(d2 * 1mm)
    sheets[5]
    substrate
    Layer() ]
steering = (θ=0, ϕ=0)
flist = 10:0.1:20

results = analyze(strata, flist, steering, logfile=devnull,
                  resultfile=devnull, showprogress=false)
nothing # hide
````

The PSSFSS run took about 55 seconds on my machine.  Here are plots of the five sheets:

````@example cpss2
using Plots
default()
ps = []
for k in 1:5
    push!(ps, plot(sheets[k], unitcell=true, title="Sheet $k", linecolor=:red))
end
plot(ps..., layout=5)
savefig("cpssb1.png"); nothing  # hide
````

![](cpssb1.png)

Notice that for all 5 sheets, the unit cell is a square of constant side length and is unrotated.
We can see from the log file (of a previous run where it was not suppressed) that this allows
PSSFSS to use additional modes in the GSM cascading procedure:

```
Starting PSSFSS 1.0.0 analysis on 2022-09-14 at 13:22:05.118
Julia Version 1.8.1
Commit afb6c60d69a (2022-09-06 15:09 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
  CPU: 8 × Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.1 (ORCJIT, skylake)
  Threads: 8 on 8 virtual cores
  BLAS: LBTConfig([ILP64] libopenblas64_.so)



Dielectric layer information...

 Layer  Width  units  epsr   tandel   mur  mtandel modes  beta1x  beta1y  beta2x  beta2y
 ----- ------------- ------- ------ ------- ------ ----- ------- ------- ------- -------
     1    0.0000  mm    1.00 0.0000    1.00 0.0000     2  1208.3    -0.0    -0.0  1208.3
     2    0.1270  mm    2.17 0.0009    1.00 0.0000     0     0.0     0.0     0.0     0.0
 ==================  Sheet   1  ========================  1208.3    -0.0    -0.0  1208.3
     3    3.8100  mm    1.04 0.0017    1.00 0.0000    10  1208.3    -0.0    -0.0  1208.3
     4    0.1270  mm    2.17 0.0009    1.00 0.0000     0     0.0     0.0     0.0     0.0
 ==================  Sheet   2  ========================  1208.3    -0.0    -0.0  1208.3
     5    2.6100  mm    1.04 0.0017    1.00 0.0000    18  1208.3    -0.0    -0.0  1208.3
 ==================  Sheet   3  ========================  1208.3    -0.0    -0.0  1208.3
     6    0.1270  mm    2.17 0.0009    1.00 0.0000     0     0.0     0.0     0.0     0.0
     7    2.6100  mm    1.04 0.0017    1.00 0.0000    18  1208.3    -0.0    -0.0  1208.3
     8    0.1270  mm    2.17 0.0009    1.00 0.0000     0     0.0     0.0     0.0     0.0
 ==================  Sheet   4  ========================  1208.3    -0.0    -0.0  1208.3
     9    3.8100  mm    1.04 0.0017    1.00 0.0000    10  1208.3    -0.0    -0.0  1208.3
 ==================  Sheet   5  ========================  1208.3    -0.0    -0.0  1208.3
    10    0.1270  mm    2.17 0.0009    1.00 0.0000     0     0.0     0.0     0.0     0.0
    11    0.0000  mm    1.00 0.0000    1.00 0.0000     2  1208.3    -0.0    -0.0  1208.3
...
```

Layers 3 and 9 were assigned 10 modes each.  Layers 5 and 7, being thinner were assigned
18 modes each. The numbers of modes are determined automatically by PSSFSS to ensure
accurate cascading.

Here are comparison plots of PSSFSS versus highly converged CST predictions digitized from
plots presented in the paper:

````@example cpss2
using Plots, DelimitedFiles
RLl = -extract_result(results, @outputs s11db(l,l))
AR11l = extract_result(results, @outputs ar11db(l))
IL21r = -extract_result(results, @outputs s21db(r,r))
AR21r = extract_result(results, @outputs ar21db(r))

default(lw=2, xlabel="Frequency (GHz)", xlim=(10,20), xtick=10:2:20,
        framestyle=:box, gridalpha=0.3)

plot(flist,RLl,title="LHCP → LHCP Return Loss", label="PSSFSS",
         ylabel="Return Loss (dB)", ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/ericsson_cpss_digitized_rllhcp.csv", ',')
plot!(cst[:,1], cst[:,2], label="CST")
savefig("cpssb2.png"); nothing  # hide
````

![](cpssb2.png)

````@example cpss2
plot(flist,AR11l,title="LHCP → LHCP Reflected Axial Ratio", label="PSSFSS",
         ylabel="Axial Ratio (dB)", ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/ericsson_cpss_digitized_arlhcp.csv", ',')
plot!(cst[:,1], cst[:,2], label="CST")
savefig("cpssb3.png"); nothing  # hide
````

![](cpssb3.png)

````@example cpss2
plot(flist,AR21r,title="RHCP → RHCP Transmitted Axial Ratio", label="PSSFSS",
     ylabel="Axial Ratio (dB)", ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/ericsson_cpss_digitized_arrhcp.csv", ',')
plot!(cst[:,1], cst[:,2], label="CST")
savefig("cpssb4.png"); nothing  # hide
````

![](cpssb4.png)

````@example cpss2
plot(flist,IL21r,title="RHCP → RHCP Insertion Loss", label="PSSFSS",
         ylabel="Insertion Loss (dB)", ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/ericsson_cpss_digitized_ilrhcp.csv", ',')
plot!(cst[:,1], cst[:,2], label="CST")
savefig("cpssb5.png"); nothing  # hide
````

![](cpssb5.png)

Differences between the PSSFSS and CST predictions are attributed to the fact that the
metalization thickness of 18 μm was included in the CST model but cannot be accommodated by PSSFSS.

