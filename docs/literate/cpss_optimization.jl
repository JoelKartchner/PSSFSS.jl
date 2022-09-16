#nb # %% A slide [markdown] {"slideshow": {"slide_type": "Slide"}}
# ## Meanderline-Based CPSS
# A "CPSS" is a circular polarization selective structure, i.e., a structure that reacts differently
# to the two senses of circular polarization.  
# We'll first look at analyzing a design presented in the literature, and then proceed to optimize another 
# design using PSSFSS as the analysis engine inside the optimization objective function.
# ### Sjöberg and Ericsson Design
# This example comes from the paper D. Sjöberg and A. Ericsson, "A multi layer meander line circular 
# polarization selective structure (MLML-CPSS)," The 8th European Conference on Antennas and Propagation 
# (EuCAP 2014), 2014, pp. 464-468, doi: 10.1109/EuCAP.2014.6901792.
#
# The authors describe an ingenious structure consisting of 5 progressively rotated meanderline sheets, which
# acts as a circular polarization selective surface: it passes LHCP (almost) without attenuation or 
# reflection, and reflects RHCP (without changing its sense!) almost without attenuation or transmission.
#

# Here is the script that analyzes their design:

using PSSFSS
## Define convenience functions for sheets:
outer(rot) = meander(a=3.97, b=3.97, w1=0.13, w2=0.13, h=2.53+0.13, units=mm, ntri=600, rot=rot)
inner(rot) = meander(a=3.97*√2, b=3.97/√2, w1=0.1, w2=0.1, h=0.14+0.1, units=mm, ntri=600, rot=rot)
center(rot) = meander(a=3.97, b=3.97, w1=0.34, w2=0.34, h=2.51+0.34, units=mm, ntri=600, rot=rot)
## Note our definition of `h` differs from that of the reference by the width of the strip.
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
##
results = analyze(strata, flist, steering, showprogress=false, 
                  resultfile=devnull, logfile=devnull); 

# Here are plots of the five meanderline sheets:

using Plots
plot(outer(rot0), unitcell=true, title="Sheet1")
#-
plot(inner(rot0-45), unitcell=true, title="Sheet2")
#-
plot(center(rot0-2*45), unitcell=true, title="Sheet3 (Center)")
#-
plot(inner(rot0-3*45), unitcell=true, title="Sheet4")
#-
plot(outer(rot0-4*45), unitcell=true, title="Sheet5")

# Notice that not only are the meanders rotated, but so too are the unit cell rectangles.
# This is because we used the generic `rot` keyword argument that rotates the entire unit
# cell and its contents. `rot` can be used for any FSS or PSS element type.  As a consequence of
# the different rotations applied to each unit cell, interactions between sheets due to higher-order
# modes cannot be accounted for; only the dominant ``m=n=0`` TE and TM modes are used in cascading
# the individual sheet scattering matrices.  This approximation is adequate for sheets that are 
# sufficiently separated.  We can see from the log file (saved from a previous run where it was not
# disabled) that only 2 modes are used to model the interactions between sheets:
#
# ```
# Starting PSSFSS 1.0.0 analysis on 2022-09-12 at 16:35:20.914
# Julia Version 1.8.1
# Commit afb6c60d69 (2022-09-06 15:09 UTC)
# Platform Info:
#   OS: Windows (x86_64-w64-mingw32)
#   CPU: 8 × Intel(R) Core(TM) i7-9700 CPU @ 3.00GHz
#   WORD_SIZE: 64
#   LIBM: libopenlibm
#   LLVM: libLLVM-13.0.1 (ORCJIT, skylake)
#   Threads: 8 on 8 virtual cores
#   BLAS: LBTConfig([ILP64] libopenblas64_.dll)
# 
# 
# 
# ******************* Warning ***********************
#    Unequal unit cells in sheets 1 and 2
#    Setting #modes in dividing layer 3 to 2
# ******************* Warning ***********************
# 
# 
# ******************* Warning ***********************
#    Unequal unit cells in sheets 2 and 3
#    Setting #modes in dividing layer 5 to 2
# ******************* Warning ***********************
# 
# 
# ******************* Warning ***********************
#    Unequal unit cells in sheets 3 and 4
#    Setting #modes in dividing layer 7 to 2
# ******************* Warning ***********************
# 
# 
# ******************* Warning ***********************
#    Unequal unit cells in sheets 4 and 5
#    Setting #modes in dividing layer 9 to 2
# ******************* Warning ***********************
# 
# 
# Dielectric layer information... 
# 
#  Layer  Width  units  epsr   tandel   mur  mtandel modes  beta1x  beta1y  beta2x  beta2y
#  ----- ------------- ------- ------ ------- ------ ----- ------- ------- ------- -------
#      1    0.0000  mm    1.00 0.0000    1.00 0.0000     2  1582.7    -0.0    -0.0  1582.7
#  ==================  Sheet   1  ========================  1582.7    -0.0    -0.0  1582.7
#      2    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#      3    4.0000  mm    1.05 0.0000    1.00 0.0000     2  1582.7    -0.0    -0.0  1582.7
#  ==================  Sheet   2  ========================   791.3  -791.3  1582.7  1582.7
#      4    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#      5    2.4500  mm    1.05 0.0000    1.00 0.0000     2   791.3  -791.3  1582.7  1582.7
#  ==================  Sheet   3  ========================     0.0 -1582.7  1582.7     0.0
#      6    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#      7    2.4500  mm    1.05 0.0000    1.00 0.0000     2     0.0 -1582.7  1582.7     0.0
#  ==================  Sheet   4  ========================  -791.3  -791.3  1582.7 -1582.7
#      8    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#      9    4.0000  mm    1.05 0.0000    1.00 0.0000     2  -791.3  -791.3  1582.7 -1582.7
#  ==================  Sheet   5  ======================== -1582.7     0.0     0.0 -1582.7
#     10    0.1000  mm    2.60 0.0000    1.00 0.0000     0     0.0     0.0     0.0     0.0
#     11    0.0000  mm    1.00 0.0000    1.00 0.0000     2 -1582.7     0.0     0.0 -1582.7
# 
# ...
# ```
# 
# Note that PSSFSS prints warnings to the log file where it is forced to set the number of layer
# modes to 2 because of unequal unit cells.  Also, in the dielectric layer list it can be seen
# that these layers are assigned 2 modes each.  The thin layers adjacent to sheets are assigned 0
# modes because these sheets are incorporated into so-called "GSM blocks" or "Gblocks" wherein
# the presence of the thin layer is accounted for using the stratified medium Green's functions.


# Here is the script that compares PSSFSS predicted performance with very
# high accuracy predictions from CST and COMSOL that were digitized from figures in the paper.

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
#md savefig("cpssa1.png"); nothing  # hide
#-
#md # ![](cpssa1.png)
#-
p = plot(flist,AR11r,title="RHCP → RHCP Reflected Axial Ratio", label="PSSFSS",
         xlim=(10,20), ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_ar_reflected.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_ar_reflected.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, ARgoal, color=:black, lw=4, label="Goal")
#md savefig("cpssa2.png"); nothing  # hide
#-
#md # ![](cpssa2.png)
#-          
p = plot(flist,IL21L,title="LHCP → LHCP Insertion Loss", label="PSSFSS",
         xlim=(10,20), ylim=(0,2), ytick=0:0.25:2)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_il.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_il.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, ILgoal, color=:black, lw=4, label="Goal")
#md savefig("cpssa3.png"); nothing  # hide
#-
#md # ![](cpssa3.png)
#-
p = plot(flist,AR21L,title="LHCP → LHCP Transmitted Axial Ratio", label="PSSFSS",
         xlim=(10,20), ylim=(0,3), ytick=0:0.5:3)
cst = readdlm("../src/assets/cpss_cst_fine_digitized_ar_transmitted.csv", ',')
plot!(p, cst[:,1], cst[:,2], label="CST")
comsol = readdlm("../src/assets/cpss_comsol_fine_digitized_ar_transmitted.csv", ',')
plot!(p, comsol[:,1], comsol[:,2], label="COMSOL")
plot!(p, foptlimits, ARgoal, color=:black, lw=4, label="Goal")
#md savefig("cpssa4.png"); nothing  # hide
#-
#md # ![](cpssa4.png)


# The PSSFSS results generally track well with the high-accuracy solutions, but are less accurate
# especially at the high end of the band, possibly because in PSSFSS metallization thickness is 
# neglected and cascading is performed for this structure using only the two principal Floquet modes.  
# As previosly discussed, this is necessary because the rotated meanderlines are achieved by rotating 
# the entire unit cell, and the unit cell for sheets 2 and 4 are not square.  Since the periodicity of 
# the sheets in the structure varies from sheet to sheet, higher order Floquet modes common to neighboring
# sheets cannot be defined, so we are forced to use only the dominant (0,0) modes which are independent of
# the periodicity.  This limitation is removed in a later example.
# Meanwhile, it is of interest to note that their high-accuracy runs
# required 10 hours for CST and 19 hours for COMSOL on large engineering workstations.  The PSSFSS run 
# took about 50 seconds on my desktop machine.


# ### Design Based on PSSFSS Optimization with CMAES
# Here we use PSSFSS in conjunction with the CMAES optimizer from the 
# [CMAEvolutionStrategy](https://github.com/jbrea/CMAEvolutionStrategy.jl) package.  I've used CMAES
# in the past with good success on some tough optimization problems.  Here is the code that defines 
# the objective function:

# ```julia
# using PSSFSS
# using Dates: now
# 
# let bestf = typemax(Float64)
#     global objective
#     """
#         result = objective(x)
# 
#     """
#     function objective(x)
#         ao, bo, ai, bi, ac, bc, wo, ho, wi, hi, wc, hc, t1, t2 = x
#         (bo > ho > 2.1*wo && bi > hi > 2.1*wi && bc > hc > 2.1*wc) || (return 5000.0)
# 
#         outer(rot) = meander(a=ao, b=bo, w1=wo, w2=wo, h=ho, units=mm, ntri=400, rot=rot)
#         inner(rot) = meander(a=ai, b=bi, w1=wi, w2=wi, h=hi, units=mm, ntri=400, rot=rot)
#         center(rot) = meander(a=ac, b=bc, w1=wc, w2=wc, h=hc, units=mm, ntri=400, rot=rot)
# 
#         substrate = Layer(width=0.1mm, epsr=2.6)
#         foam(w) = Layer(width=w, epsr=1.05)
#         rot0 = 0
# 
#         strata = [
#                 Layer()
#                 outer(rot0)
#                 substrate
#                 foam(t1*1mm)
#                 inner(rot0 - 45)
#                 substrate
#                 foam(t2*1mm)
#                 center(rot0 - 2*45)
#                 substrate
#                 foam(t2*1mm)
#                 inner(rot0 - 3*45)
#                 substrate
#                 foam(t1*1mm)
#                 outer(rot0 - 4*45)
#                 substrate
#                 Layer() ]
#         steering = (θ=0, ϕ=0)
#         flist = 11.5:0.25:18.5
#         resultfile = logfile = devnull
#         showprogress = false
#         results = analyze(strata, flist, steering; showprogress, resultfile, logfile)
#         s11rr, s21ll, ar11db, ar21db = eachcol(extract_result(results, 
#                        @outputs s11db(R,R) s21db(L,L) ar11db(R) ar21db(L)))
#         RL = -s11rr
#         IL = -s21ll
#         RLgoal, ILgoal, ARgoal  = (0.4, 0.5, 0.6)
#         obj = max(maximum(RL) - RLgoal, 
#                   maximum(IL) - ILgoal, 
#                   maximum(ar11db) - ARgoal, 
#                   maximum(ar21db) - ARgoal)
#         if obj < bestf
#             bestf = obj
#             open("optimization_best.log", "a") do fid
#                 xround = map(t  -> round(t, digits=4), x)
#                 println(fid, round(obj,digits=5), " at x = ", xround, "  #", now())
#             end
#         end
#         return obj
#     end
# end
# ```

# We optimize at 29 frequencies between 11.5 and 18.5 GHz.  In the previously presented 
# Sjöberg and Ericsson design a smaller frequency range of 12 to 18 GHz was used for optimization.
# We have also adopted more ambitious goals for return loss, insertion loss, and axial ratio
# of 0.4 dB, 0.5 dB, and 0.6 dB, respectively. These should be feasible because for
# our optimization setup, we have relaxed the restriction that all unit cells must be
# identical squares.  With more "knobs to turn" (i.e., a larger search space), we expect to be able
# to do somewhat better in terms of bandwidth and worst-case performance.  For our setup there are 
# fourteen optimization variables in all.  
# As one can see from the code above, each successive sheet in the structure is rotated an additional
# 45 degrees relative to its predecessor.   The objective is defined as the maximum departure
# from the goal of RHCP reflected return loss, LHCP insertion loss, or reflected or transmitted axial ratio that
# occurs at any of the analysis frequencies (i.e. we are setting up for "minimax" optimization). Also,
# the `let` block allows the objective function to maintain persistent state in the
# variable `bestf` which is initialized to the largest 64-bit floating point value. Each time a set 
# of inputs results in a lower objective function value, `bestf` is updated with this value and 
# the inputs and objective function value are
# written to the file "optimization_best.log", along with a date/time stamp.  This allows the user
# to monitor the optimization and to terminate the it prematurely, if desired, without losing the 
# best result achieved so far. Each objective function evaluation takes about 9 seconds on my machine.

# Here is the code for running the optimization:

# ```julia
# using CMAEvolutionStrategy
# #  x = [a1,  b1,   a2, b2,  a3,  b3,  wo,  ho,  wi,   hi,  wc,   hc,  t1,  t2]
# xmin = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.1, 0.1, 0.1,  0.1, 0.1,  0.1, 2.0, 2.0]
# xmax = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 0.35,4.0, 0.35, 4.0, 0.35, 4.0, 6.0, 6.0]
# x0 = 0.5 * (xmin + xmax)
# isfile("optimization_best.log") && rm("optimization_best.log")
# popsize = 2*(4 + floor(Int, 3*log(length(x0))))
# opt = minimize(objective, x0, 1.0;
#            lower = xmin,
#            upper = xmax,
#            maxfevals = 9000,
#            popsize = popsize,
#            xtol = 1e-4,
#            ftol = 1e-6)
# ```

# Note that I set the population size to twice the normal default value.  Based
# on previous experience, using 2 to 3 times the default population size helps the 
# optimizer to do better on tough objective functions like the present one.
# I let the optimizer run overnight for about 17 hours, during which time it reduced the objective function
# value from 11.48 dB to -0.12 dB. Since it appeared to have essentially converged, I terminated the run.

# Here is a look at the beginning and final portions of the file "optimization_best.log":
# ```
# 11.47502 at x = [4.3371, 3.6377, 4.6978, 3.0065, 4.2178, 3.7151, 0.3456, 1.4562, 0.3345, 2.1578, 0.3498, 2.0437, 3.3593, 3.5733]  #2022-09-13T14:35:39.320
# 8.83597 at x = [4.3381, 4.8048, 3.4625, 4.7309, 3.051, 3.5356, 0.1072, 3.585, 0.1039, 2.2732, 0.3206, 1.1203, 3.8363, 2.2236]  #2022-09-13T14:35:57.550
# 6.61653 at x = [3.4206, 3.1036, 3.0215, 3.2064, 3.0201, 3.4882, 0.1034, 1.5586, 0.103, 2.3935, 0.2026, 2.1608, 4.2708, 2.4485]  #2022-09-13T14:36:06.768
# 6.17284 at x = [4.6927, 3.9841, 4.6078, 3.8153, 3.2699, 3.0023, 0.3254, 3.0281, 0.3167, 1.0311, 0.3334, 1.3243, 4.9482, 2.5076]  #2022-09-13T14:37:01.227
# 4.14975 at x = [4.4617, 4.105, 3.9944, 4.8857, 4.6484, 3.6307, 0.3457, 2.0488, 0.3354, 1.9338, 0.3275, 1.5484, 3.2037, 2.683]  #2022-09-13T14:37:17.647
# 2.2661 at x = [4.7303, 4.1483, 3.4428, 3.3463, 3.607, 3.2735, 0.1572, 2.2432, 0.35, 2.2474, 0.2468, 1.376, 2.6822, 2.3709]  #2022-09-13T14:39:58.103
# 1.57851 at x = [4.95, 4.1857, 4.0174, 3.3486, 4.9004, 3.678, 0.12, 1.6691, 0.2451, 1.0522, 0.2964, 1.7673, 3.6622, 2.6612]  #2022-09-13T14:41:33.641
# ...
# -0.11824 at x = [3.0642, 4.8748, 4.1393, 3.0066, 3.8585, 3.0047, 0.344, 3.0002, 0.3469, 1.4042, 0.2502, 2.3751, 3.766, 2.3372]  #2022-09-14T06:37:45.783
# -0.1184 at x = [3.0416, 4.8731, 4.148, 3.0012, 3.8553, 3.0103, 0.3449, 3.0104, 0.3482, 1.4088, 0.2533, 2.3858, 3.7593, 2.3347]  #2022-09-14T06:38:57.587
# -0.11896 at x = [3.0454, 4.881, 4.1304, 3.0001, 3.8389, 3.01, 0.345, 3.0042, 0.348, 1.413, 0.2524, 2.3719, 3.7682, 2.3398]  #2022-09-14T06:50:21.349
# -0.11961 at x = [3.069, 4.8773, 4.1652, 3.002, 3.8446, 3.0055, 0.3439, 3.0101, 0.3468, 1.4051, 0.2518, 2.3774, 3.7652, 2.3401]  #2022-09-14T06:59:09.133
# -0.11961 at x = [3.0661, 4.8788, 4.155, 3.0014, 3.8314, 3.0069, 0.3441, 2.9933, 0.3476, 1.4168, 0.2516, 2.3689, 3.7816, 2.3322]  #2022-09-14T07:06:13.332
# -0.1197 at x = [3.0657, 4.8757, 4.1549, 3.0018, 3.8331, 3.0061, 0.3436, 3.0024, 0.3468, 1.4055, 0.2517, 2.3716, 3.7758, 2.3382]  #2022-09-14T07:15:07.056
# ```

# The final sheet geometries and performance of this design are shown below:

# ![](./assets/cpss_cmaesopt_sheets.png)

# ![](./assets/cpss_cmaesopt_rl_refl.png)

# ![](./assets/cpss_cmaesopt_ar_refl.png)

# ![](./assets/cpss_cmaesopt_il_trans.png)

# ![](./assets/cpss_cmaesopt_ar_trans.png)

# As hoped for, the performance meets the more stringent design goals over a broader bandwidth than the 
# Sjöberg and Ericsson design, presumably because of the greater design flexibility allowed here.
