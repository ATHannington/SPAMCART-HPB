# SPAMCART-HPB
### **SPAMCART-HPB** - An independent, partial implementation of the SPAMCART method, which combines Smoothed Particle Hydrodynamics and Monte Carlo Radiative Transfer. 

SPAMCART addresses the problem of applying Monte Carlo Radiative Transfer (MCRT) directly to Smoothed Particle Hydrodynamics (SPH) without an intervening step wherein the SPH fluid (where fluid is modelled with large numbers of particles of a fixed mass) is mapped to 3D grid. Mapping to such a grid severely limits the final simulation resolution, especially if MCRT is performed "on-the-fly" as the hydrodynamical (SPH) simulation is run (rather than the MCRT applied solely as a post-processing step). In such "on-the-fly" radiative-hydrodynamics simulations, the resolution is reduced when converting from SPH to a grid for MCRT processing **_and_** is reduced again when converting the MCRT grid back to SPH particles for the next simulation timestep.

Developed as a methodological validation project during my undergraduate studies at the University of St Andrews, SPAMCART-HPB is a partial implementation of the SPAMCART method. The software presented here was developed independently from the original SPAMCART team, to independently derive and verify the original [O. Lomax & A. P. Whitworth 2016, SPAMCART method](http://cdsads.u-strasbg.fr/abs/2016MNRAS.461.3542L).

## Affiliations, credits, and acknowledgements

*Author*: Andrew Tomos Hannington
    
*Affiliation*: University of St Andrews

*Credit*: SPAMCART-HPB is a direct derivative of the work by [O. Lomax & A. P. Whitworth 2016](http://cdsads.u-strasbg.fr/abs/2016MNRAS.461.3542L) [SPAMCART algorithm](https://github.com/odlomax/spamcart-dev). The software presented in this repository is an independent, partial implementation of the original Lomx and Whitworth SPAMCART algorithm.

*Acknowledgements*: Work by Andrew Hannington as part of Master's project work with Dr Maya Petkova, and Professor Ian Bonnell at the University of St Andrews.

