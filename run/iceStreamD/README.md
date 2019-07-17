# Simulations to investigate the spatial development of drainage on Bindschader Ice Stream
This directory contains the simulations for Section 4.1 and Figure 2 in *Elsworth and Suckale* (2016) to determine the evolution of basal strength along the shear margin in Bindschadler Ice Stream. In order to match observed surface velocity data, we find that a concentration of basal strength is necessary in the shear margin. Assuming the presence of efficient drainage in this region, we can connect these best fitting basal strength profiles to the relative contribution of efficient and distributed meltwater drainage. We find that along flow increased discharge occurs in an efficient drainage element underneath the margin, while the distributed system becomes less well-connected. 

### Contents
*plotHydroProperties.m*: Plot the discharge and film thickness of an assumed Rothlisberger channel, thin film hydrological system. This script reproduces Figures 2d and 2e from *Elsworth and Suckale* (2016).

*plotStrength.m*: Plot the best fitting strength distributions to match surface velocity data. This script reproduces Figure 2b from *Elsworth and Suckale* (2016).

*plotVelocity.m*: Plot the best-fitting velocity model against the observed surface velocity along Bindschadler Ice Stream. This script reproduces Figure 2a from *Elsworth and Suckale* (2016).

*profile_/* : Simulation set-up and results for the profiles along Bindschadler Ice Stream. The locations of these profiles are shown in Figure 1b of *Elsworth and Suckale* (2016).

*profileVelData.mat*: Observational surface velocity data from *Scambos et al.* (1994) and *Bindschadler et al.* (1996). The locations of these profiles are shown in Figure 1b and velocities are shown in Figure 2a [*Elsworth and Suckale, 2016*].

### References
Bindschadler R A, P Vornberger, D D Blankenship, T A Scambos, and R W Jacobel (1996) Surface velocity and mass balance of Ice Streams D and E, West Antarctica. Journal of Glaciology 42(142):461–475.

Elsworth C W and J Suckale (2016) Rapid ice flow rearrangement induced by subglacial drainage in West Antarctica. Geophysical Research Letters 43.

Scambos T A, K A Echelmeyer, M A Fahnestock, and R A Bindschadler (1994) Development of enhanced ice flow at the southern margin of Ice Stream D, Antarctica. Annals of Glaciology 20(1):313–318.
