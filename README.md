# kgh_angcorr_analysis_base
The state of the particle decay angular correlation analysis code as of my (K. Hanselman's) graduation from FSU, ca. 2022

***Edited and committed to github in August 2025, for the use of academic collaborators.

DISCLAIMER: THIS IS NOT INTENDED to be a "complete" product. It is the result of my graduate school research work, as documented in my 2022 FSU dissertation and upcoming PRC. These are all different chunks of code gradually tacked on to each other: some my own, some given by Ingo Wiedenhoever, my advisor, from his Cologne days, some the generalized C++ versions of the original fortran code given to me by Prof. Jeffrey Tostevin of U. Surrey. What I've done here is tried to clean it up and pare it down as much as possible so that it can be downloaded, manipulated, and improved by any collaborators interested in taking it farther.

What you will need (prereq's):
  >> C++ compiler
  >> ROOT
  >> FRESCO
  >> GSL

The last time I compiled and ran the codes was as of this upload, on the Ubuntu partition of my Windows 11 personal laptop. The ROOT version was 6.28/04, g++ 11.4.0, and GSL 2.7.1

STEP 1: FRESCO

First, a reaction calculation must be performed in FRESCO, populating the excited state which will decay. This can be any reaction, and with either gamma or nucleon decays in mind, but nothing heavier. The math as implemented is only for gammas, protons, or neutrons. The parameter LAMPL on card 5 should be set to -2, to print out the scattering amplitudes on fort.37 for the exit channel of a typical transfer reaction. This file needs to be saved in the "37files/" directory under some name. *There needs to be ONE excited state per file; the code from here on will only take the first state, so multi-state outputs need to be split up manually.

STEP 2: Decay Coupling

Next, the density matrices and orientation tensors are calculated from the fort.37 files, and the appropriate decay coefficients are coupled to them with the spherical harmonics to form the decay patterns. Which code you use for this depends on the physical situation.

All take text file inputs as: ./exec input_file.txt, which have their own directory.

>> bees.cpp -- nucleon decay from a single unbound state
 
   In this case, inputs are taken from "bees_inputs/", where the relevant parameters
   of the reaction are defined. This includes the fort.37 file to use, the parity
   and branching ratio, and the Euler angles used to rotate the coordinate system
   from the default Madison coordinates of FRESCO, where z is along the beam axis.
   Many things are then calculated: the 2D decay patterns, the m-substate distributions
   (along the beam axis), the density matrices and orientation tensors themselves, etc.
   Because the coordinate system is variable, the angle definitions are made general:
   "xy" angle is that in the x-y plane, e.g. "phi", while the "z" angle is that relative
   to the z-axis, e.g. "theta" -- typically in cosine. These calculations are dependent
   on many external codes: "WikoParticleClass.h" for the decay coefficients, "BeesHolder.h"
   for the matrix/tensor formation, and "wiko_wiko5.hpp" for the mathematical functions.

   The idea is to keep this step as simple and general as possible so that its outputs
   can then be used for more complicated analysis. Keep things "modular".

>> gambees.cpp -- gamma decay from a bound state

   This proceeds nearly identically to bees.cpp, except that it's gamma decay from a
   bound state. The only reason it deserves its own code is because the coefficients
   are different, hence requiring a different code, "WikoGammaClass.h". It seemed
   cleaner to do things this way.

>> geebees.cpp -- nucleon decay from multiple OVERLAPPING unbound states

   The reason this is different from bees.cpp is because it requires a "grand" density
   matrix due to the multiple initial states. Hence it needs only one code, "GrandBeesHolder.h",
   which has everything all included in it (because it's all coupled so early on).

   It requires MULTIPLE fort.37 files as input, one for each initial state. These must each
   be in their OWN 37 file as explained earlier. Besides the usual parity and branching ratio,
   we also need to know the excitation energy and width as these potentially enter in as extra
   factors. See my dissertation, upcoming PRC, or (if you're really desperate) my research logs
   at FSU for a more detailed discussion of this physics.

   Apart from keeping track of "individual" versus "coupled" parts of the calculation, however,
   this proceeds pretty much like bees.cpp, with its outputs intended to be used later
   for other codes.


Rootfile outputs are stored appropriately in "rootfiles". These rootfiles are what are taken in by later codes
as inputs.

As far as these "later codes" are concerned, I've left in two examples of what I used to do:

 >> phi_decomp_analysis: this is where the "phi" (xy) dimension, typically taken in Madison coordinates, is
                         transformed into a discrete "q" index, turning the full 2D distribution into multiple
                         1D distributions. This analysis was motivated by the 2020 experiments with SABRE which
                         has a full phi dimension but incomplete costheta. The inputs are all fixed for this,
                         as the maximum q to transform is fixed by the spin: qmax = kmax = 2J. The main outputs
                         for this are text files which can be easily plotted.

 >> projection_analysis: this analysis is poorly named, because it's more 1D "slices" of the distributions than
                         it is a true projection. Therefore the slices themselves actually do need to be defined,
                         similar to the 2021 barrel data only without bin widths. Then comparison can be made to
                         the data in a similar way as the decomps.


Again, this is not a finished product, just the state of things when I laid this kind of analysis down. Any collaborators wishing to pick it up are more than welcome to change it, improve it, do things differently, etc. You can reach me at kghanselman@gmail.com if there are any questions.

-- K. Hanselman | August 2025
