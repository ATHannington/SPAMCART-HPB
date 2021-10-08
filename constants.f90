!-------------------------------------------------------------------------------------
module constants
!=============================================
!	MCRT Constants:
!=============================================
		integer,save :: Nmc=1000000!20000000!200000				!Number of monte carlo packets

		integer :: seed=1234							!Seed for random numbers
		integer :: nBins=500							!Number of bins for bin arrays

		character(len=100) ::filename = "tau_scatters_136_large_pt2.dat"!"mcrt_136.dat"		!Output filename

		Real*8,allocatable,dimension(:,:) :: bin				!Bin for Next Event Estimator Data output
		Real*8 :: binWidth							!Width of bins, set in MCRT main
		Real*8,dimension(3) :: r_im = (/0.d0,0.d0,2.d0/)			!Position vector of image plane
		Real*8 :: theta_im,phi_im						!Angles of image plane
		Real*8 :: systemwidth = 2.0d0						!Width of square box to bin over

		Real*8 :: vacuum_tol = 1d-5						!Vacuum Tol sets additional range that s0 can take to not be in vacuum, to account for precision errors.
		
		Integer :: inverseprintpercentage = 1000				!Sets inverse of percentage to print completion for mc packets. i.e. =100 => every 1% will be printed

		integer :: loopmax = 10000						!Sets maximum number of iterations loops can reach before an infinite loop is detected
											!Set this very high or comment out relevant code in MCRT if high optical depths being used.
	!=============================================
	!	MCRT multiple sources Constants:
	!=============================================
		integer :: n_sources							!Number of radiation sources
		Real*8 :: Lumi = 1.d0							!Luminosity of current radiation source (power) [W]	
		Real*8, allocatable,dimension(:) :: Lumi_array				!Array of source Luminosities for all sources(power) [W]
		Real*8, allocatable, dimension(:,:) :: source_origin_array		!Array containing the origins of all sources
		Real*8,dimension(3) :: source_origin = (/0.d0,0.0d0,0.d0/)		!Updated in MCRT to look at current source location, but can also be used as fixed source location if
											!n_sources =1 and "call initialise_sources" is commented out in MCRT main.

		character(len=50) ::radsourcefile = "Rad-Sources-xyzL-2-Even-Sources.dat"	!Filename for radiation source data. Format is x,y,z, luminosity with a header line of n_sources

!=============================================
!	SPH Constants:
!=============================================

		character(len=40) ::sphfile = "uniform-sphere-xyzhm-136.dat"		!SPH particle data file

		Real*8 :: chi = 10.d0							!Mass specific absorption coefficient
		Real*8 :: tauchim1							!tau*(chi**(-1))

		Real*8,allocatable,dimension(:) :: h_array, m_array			!Arrays for particle masses and smooting lengths

!Constant h and m for testing purposes:
		Real*8 :: hconst = 1.d0							!Constant smoothing length for debugging
		Real*8 :: mconst = 1.d0							!Constant particle mass for debugging

!Temporary h and m that will update with each particle:
		real*8 :: zeta = 2.d0							!***CONSTANT*** Traditionally zeta=2*h (as per M4 cubic Spline kernel) but this code is predominantly
											!in units of h. (L,t and origins NOT in units of h)
		Real*8 :: h,m								!Mass and smoothing length of current particle. This is updated throughout the code whenever they are used to
											!be the mass and smoothing length of the current SPH particle being used.

!=============================================
!	Vacuum Algorithm Constants:
!=============================================
		Real*8 :: impact_fraction = 5.0d-2					!The depth a photon will penetrate an SPH particle in units of h.

!=============================================
!	Dist_Finder Constants:
!=============================================
		Real*8 :: taudum							!Dummy optical depth. This varies in the distance finder only.

		Real*8 :: overest = 1.2d0!!!1.2d0					!Over-estimation of l0, parameter

		Real*8 :: tol=1.d-6							!Tolerance of distance finder. Gives how close last two guesses (lnp1 and ln) should be to 
											!be considered converged
		Real*8 :: lneg = 1.d-3							!Sets the value ln should be set to if ln has gone negative.

		Real*8 :: convergence_factor = 5.d-1					!Factor to decrease l0 by if it's too large. 
		integer :: convergeattmeptlim = 200					!Number of attempts distance finder should make to attempt convergence. Will take average of last 10
											!guesses if convergence fails
		Real*8 :: l0_min = 1d-8							!Minimum l0 guess before distance finder assumes no positive solution. Make this very small for high optical
											!depth systems!!!
						

			!***++++***
!***DO NOT SET THESE TO FALSE!!! THEY ARE FOR SINGLE PARTICLE TESTS ONLY!!! :***
		Logical :: UnderEstBool = .true.					!Turns under estimate of l0 system on
		Logical :: lnm1SideBool = .true.					!Turns on alternating side of lnm1 to ensure bracketing of root on
			!***++++***

!=============================================
!	Intercepts Constants:
!=============================================
		Real*8,allocatable,dimension (:,:) :: xi_intersects			!Holds the origins of the SPH particles the photon packet ray will intersect
		Integer :: nintersects							!Holds the number of SPH particles the photon packet ray will intersect
		Real*8,allocatable,dimension(:) :: h_intercept_array, m_intercept_array		!Arrays for particle masses and smooting lengths of the SPH particles the ray intersects
!=============================================
!	Physical Constants:
!=============================================
		Real*8 :: pi = ACOS(-1.d0)						!PI
		Real*8 :: c = 2.998d+8 							![m s^-1] Speed of Light 
		Real*8 :: AU = 1.496d+11 						![m] Astronomical Unit in meters

end module constants



!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!	Module for temporary global variables pertaining to the 
!	coordinates of current ray and SPH particle.
!-------------------------------------------------------------------------------------
module coordinates
	Real*8 :: t, c_dist, s0, sl						
!t==Closest approach **[code units]**
!c_dist==Impact parameter [code units/h] -> [dimensionless]
!s0==radial distance to start position/origin [code units/h] -> [dimensionless]
!sl== radial distance to packet end point [code units/h] -> [dimensionless]

!^^^ All of above w.r.t. SPH particle centre.

!	Temp Coordinate constants (vectors here):
	Real*8,dimension(3) :: xivec,ovec,nvec						!Temporary vectors of xi (particle origin), o(packet origin) and n-hat (packet direction)
											!Used in distance finder where some of these (especially ovec) is needed to be altered.

!	True Coordinate constants:
	Real*8,dimension(3) :: ntruevec = (/0.0d0,0.0d0,0.0d0/)				!"True" n-hat vector of packet direction as cast by the MCRT main program.
	Real*8,dimension(3) :: otruevec = (/0.0d0,0.0d0,0.0d0/)				!"True" o vector of packet origin as cast by the MCRT main program.
	Real*8,dimension(3) :: xitruevec = (/0.0d0,0.0d0,0.0d0/)			!"True" xi vector of SPH particle origin.

	!Particle Locations
	Real*8,allocatable,dimension (:,:) :: xi_array					!Array of "True" xi_vectors.
	Integer :: n_particles = 3							!This should match the second dimension of xi_array. Set in "call initialise_particles"

	Real*8,dimension(3) :: xi1 = (/0.0d0,0.d0,0.d0/)				!Dummy SPH particle positions for test purposes
	Real*8,dimension(3) :: xi2 = (/0.0d0,-2.0d0,0.0d0/)
	Real*8,dimension(3) :: xi3 = (/0.0d0,2.0d0,0.d0/)

end module coordinates
!-------------------------------------------------------------------------------------
