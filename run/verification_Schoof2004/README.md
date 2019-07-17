# Verification of the code developed in *Elsworth and Suckale* (2016) against the analytical solution of *Schoof* (2004)
We compare our numerical solution of the margin tracking algorithm to the analytical solution developed by *Schoof* (2004). The analytical solution assumes a Newtonian ice rheology and solves for the margin location consistent with the force balance. Figures 3(a1), 3(a2) and 3(a3) in *Schoof* (2004) show the assumed basal strength profiles, determined shear margin location and basal stress. The algorithm developed in *Elsworth and Suckale* (2016) recovers the same margin position and basal strength profile as the analytical solution and is therefore verified for a Newtonian rheology. The inclusion of a non-Newtonian rheology does not change the solution technique for determining the shear margin location, so we posit that the algorithm remains accurate.

# Contents
*data/*: Digitized analytical solution profiles from Figures 3(a1), 3(a2) and 3(a3) in *Schoof* (2004).
*figure3_1/*: Simulation set-up and results for the basal strength profile in Figure 3(a1) of *Schoof* (2004).  
*figure3_2/*: Simulation set-up and results for the basal strength profile in Figure 3(a2) of *Schoof* (2004). This simulation is included in *Elsworth and Suckale* (2016) in Figure S3.
*figure3_3/*: Simulation set-up and results for the basal strength profile in Figure 3(a3) of *Schoof* (2004).
*plotBasalStressVerification*: Function to plot the model result and analytical solution. 
*plotSchoofVerification*: Plots the model simulations against the analytical solution. Reproduces *benchmrk_Schoof2004.eps* and Figure S3 from *Elsworth and Suckale* (2016).
*verification_Schoof2004.eps*: Figure comparing the model simulations and analytical solutions. 

# References
Elsworth C W and J Suckale (2016) Subglacial drainage may induce rapid ice flow rearrangement in West Antarctica. Geophysical Research Letters (in press).

Schoof C G (2004) On the mechanics of ice-stream shear margins. Journal of Glaciology 50(169):208–218.
