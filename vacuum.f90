!-----------------------------------------------------------------------------------
!Subroutine to check whether the current origin lies within any of the SPH particles
!-----------------------------------------------------------------------------------
subroutine vacuum_check(ovecdum,VacuumBool)
	use constants
	use coordinates
	implicit none
		Real*8,dimension(3),intent(in) :: ovecdum
		Logical, intent(out) :: VacuumBool
		integer :: i,j
		Real*8,dimension(3) :: xi
		Real*8 :: s0dum
		Real*8 :: dum

	VacuumBool = .TRUE.

	Do i=1,nintersects
		xi(1:3) = xi_intersects(1:3,i)

!important that we set new h and m for each xi:
		call set_particle_params(i)

		dum = 0.d0
		Do j=1,3
			dum = dum + (ovecdum(j)-xi(j))**2
		enddo
		dum = DSQRT(dum)/h
		s0dum = dum

	!If s0<=zeta then packet lies within a particle and so is NOT in a vacuum:
	!+tol to allow for rounding precision errors: 
		If(s0dum.le.(zeta+vacuum_tol)) then
			VacuumBool = .FALSE.
			GOTO 33
		endif

	enddo

	33 continue

	RETURN
end subroutine vacuum_check

!-----------------------------------------------------------------------------------
!Subroutine to give distance to nearest intercepting particle when packet is in Vacuum
!-----------------------------------------------------------------------------------
subroutine Vacuum_alg(ovecdum,L,ExitBool)
	use constants
	use coordinates
	implicit none
		Real*8,dimension(3),intent(in) :: ovecdum
		Real*8,intent(out) :: L
		Logical, intent(out) :: ExitBool
		integer :: i,j
		Real*8,dimension(3) :: xi
		Real*8,allocatable,dimension(:) :: l_array
		Real*8 :: ldum
		Logical :: SolBool
		integer :: minlocation
		Real*8 :: tdum, cdum
		Real*8 :: t_error
		Real*8 :: sldum

	allocate(l_array(n_particles))

	ldum = 1.d125
	l_array = 1.d125

!Set desired sl to be zeta-impact_fraction
	sldum = (zeta-impact_fraction)

	j = 0
	Do i=1,n_particles
		xi(1:3) = xi_array(1:3,i)
	
	!Needs masses and smoothing lengths for all particles, not just intercepts:	
		h = h_array(i)
		m = m_array(i)
!Print*,"h=",h
!Print*,"m=",m

!Print*,"[@Vac alg:] zeta",zeta
!Print*,"[@Vac alg:] impact_fraction",impact_fraction
!Print*,"[@Vac alg:] sldum",sldum

	!Get t and c, not affected by zeta:
		call get_t(xi,ovecdum,tdum)
		call get_c(tdum,xi,ovecdum,cdum)

!Print*,"[@Vac alg:] tdum",tdum
!Print*,"[@Vac alg:] cdum",cdum

	!Sets all quantities apart from sl, which is not needed here
		call get_coords(xi,ovecdum,0.d0)

	!Get distance between point of closest approach and particle edge NOT in units of h:
		t_error = DSQRT(((h*zeta)**2)-((h*c_dist)**2))

	!Sanity check:
		If(isnan(t_error).eqv..TRUE.)then
			Print*,"[@Vac alg:] t_error NaN! Program will Abort!"
			Print*,"[@Vac alg:] zeta=",zeta
			Print*,"[@Vac alg:] cdum=",cdum
			Print*,"[@Vac alg:] c_dist=",c_dist
			Print*,"[@Vac alg:] t_error=",t_error
			Print*,"[@Vac alg:] tdum=",tdum
			Print*,"[@Vac alg:] t+t_error=",t+t_error
			STOP
		endif


	!Include tdum check to ensure forward particles only
	!tdum + t_error checks for packet lying in particles behind it, but that packet is still inside.
	!C<=zeta means packet will intersect:
		If (((tdum+t_error).gt.(0.d0)).and.(cdum<=zeta)) then
		!find distance to move through "vacuum" to intercept particle:
			call get_l(sldum,xi,ovecdum,-1.d0,ldum,SolBool)
!Print*,"[@Vac alg:] ldum",ldum
		!If no solution to distance, set l to infinity to rule out as minimum distance.
		!If solution exists, add this particle to intercepts arrays:
			If(SolBool .eqv. .FALSE.) then
				ldum = 1.d125
			else if(SolBool .eqv. .TRUE.) then
				j = j + 1
				xi_intersects(1:3,j) = xi(1:3)

				h_intercept_array(j) = h
				m_intercept_array(j) = m
			endif

		!Ignore particles behind current location
			If(ldum.le.(0.d0))then
				ldum=1.d125
			endif
		endif

!Print*,"[@Vac alg:] ldum",ldum

		l_array(i) = ldum

	enddo


!number of intersected particles is equal to the counter j
	nintersects = j

!Choose nearest particle to move to:
	minlocation=minloc(l_array,1)
	L = l_array(minlocation)	

!If no particles intersected, kill the packet:
	If (nintersects.gt.(0))then
		ExitBool = .false.
	else
		ExitBool = .true.
	endif

	If(ExitBool .eqv. .TRUE.) then
		L = 0.d0
	endif

!Print*,"[@Vac alg:] nintersects",nintersects
!Print*,"[@Vac alg:] xi_array(1:3,1:5)",xi_array(1:3,1:5)
!Print*,"[@Vac alg:] xi_intersects(1:3,1:5)",xi_intersects(1:3,1:5)
!Print*,"[@Vac alg:] xi_intersects(1:3,1:nintersects)",xi_intersects(1:3,1:nintersects)
!Print*,"[@Vac alg:] ExitBool",ExitBool
!STOP
	RETURN
end subroutine Vacuum_alg

!-----------------------------------------------------------------------------------
!Subroutine to give distance L for given ovec, xivec, nvec, and sl
!-----------------------------------------------------------------------------------
subroutine get_l(sldum,xi,ovecdum,plusminus,l,SolBool)
	use constants
	use coordinates
	implicit none
		Real*8,intent(in) :: sldum
		Real*8,dimension(3),intent(in) :: xi
		Real*8,dimension(3),intent(in) :: ovecdum
		Real*8,intent(in) :: plusminus
		Real*8,intent(out) :: l
		Logical, intent(out) :: SolBool
		Real*8 :: adum,bdum,cdum
		Real*8 :: dum, mag
		integer :: i
		Real*8,dimension(3) :: svecdum
		Real*8 :: sqrt_arg

	svecdum = (ovecdum - xi)

	adum = 1.d0
!Print*,"[@get_l:] a=",adum

	dum=0.d0
	Do i=1,3
		dum = dum + (ntruevec(i)*svecdum(i))
	enddo
	bdum = dum*2.d0

!Print*,"[@get_l:] b=",bdum

	dum=0.d0
	Do i=1,3
		dum = dum + (svecdum(i))**2
	enddo
	cdum = dum - ((sldum*h)**2)
!Print*,"[@get_l:] h=",h
!Print*,"[@get_l:] plusminus=",plusminus
!Print*,"[@get_l:] sldum=",sldum
!Print*,"[@get_l:] c=",cdum
!Print*,"[@get_l:] bdum**2- 4.d0*adum*cdum=",bdum**2- 4.d0*adum*cdum
!Print*,"[@get_l:] (2.d0*adum)=",(2.d0*adum)

	sqrt_arg = bdum**2- 4.d0*adum*cdum

!Prevent issues with no solution.
!If SQRT is NaN no physical solution exists, and as such solbool === false.
	If(sqrt_arg  .lt. (0.d0)) then
		SolBool = .FALSE.
	else
		l = ((-1.d0*bdum) + (plusminus*DSQRT(sqrt_arg))) / (2.d0*adum)
	endif

!Print*,"[@get_l:] l=",l	
	If(isnan(l) .eqv. .true.) then
		SolBool = .FALSE.
	else
		SolBool = .true.
	endif

	RETURN
end subroutine get_l

!-----------------------------------------------------------------------------------
!Subroutine to give t
!-----------------------------------------------------------------------------------
subroutine get_t(xi,ovecdum,tdum)
	use constants
	use coordinates
	implicit none
		Real*8,dimension(3),intent(in) :: xi
		Real*8,dimension(3),intent(in) :: ovecdum
		Real*8,intent(out) :: tdum
		integer :: i
		Real*8 :: dum

! Calculate interaction parameters:

	!Closest approach, t:
	i=1
	dum = 0.d0
	Do i=1,3
		dum = dum + (xi(i) - ovecdum(i))*ntruevec(i)
	enddo

	tdum = dum

	RETURN
end subroutine get_t

!-----------------------------------------------------------------------------------
!Subroutine to give c for a given t, NOT AFFECTED BY ZETA!!!
!-----------------------------------------------------------------------------------
subroutine get_c(tdum,xi,ovecdum,cdum)
	use constants
	use coordinates
	implicit none
		Real*8,dimension(3),intent(in) :: xi
		Real*8,intent(in) :: tdum
		Real*8,dimension(3),intent(in) :: ovecdum
		Real*8,intent(out) :: cdum
		integer :: i
		Real*8 :: dum

! Calculate interaction parameters:

	!Closest approach, t:
	i=1
	dum = 0.d0
	Do i=1,3
		dum = dum + (ovecdum(i) + tdum*ntruevec(i) - xi(i))**2
	enddo
	
	dum = DSQRT(dum)/h

	cdum = dum

	RETURN
end subroutine get_c
