!-----------------------------------------------------------------------
!		Subroutine of SPH particle density
!	Takes contributions from all intersecting SPH particles
!	and uses non min(s,zeta) form to ensure no errors with s0=zeta:
!-----------------------------------------------------------------------
subroutine density (vec,rho)
	use constants
	use coordinates
		real*8,dimension(3),intent(in) :: vec
		Real*8,intent(out) :: rho
		real*8,dimension(3) :: xi
		integer :: i,j
		Real*8 :: dum
		Real*8 :: w, wdum
		Real*8 :: s
		Logical :: Logicaldum
		integer :: intersect_index

	wdum = 0.d0
	w = 0.d0

!Print*,"[@density:] ovec=",ovec(1:3)

	Do i=1,nintersects
		xi(1:3) = xi_intersects(1:3,i)

	!important that we set new h and m for each xi:
		call set_particle_params(i)

		dum = 0.d0
		Do j=1,3
			dum = dum + (vec(j) - xi(j))**2
		enddo
		s = DSQRT(dum)/h

		call kernel(s,wdum)
		w = w + (m/h**3)*wdum

!Print*,"[@density:] wdum=",wdum			
	enddo
	rho = w

!Print*,"[@density:] rho=",rho

	RETURN
end subroutine density

!-----------------------------------------------------------------------
!	Subroutine to normalise unit vectors
!-----------------------------------------------------------------------
subroutine normalise
	use constants
	use coordinates
	implicit none
		Real*8 :: mag
		integer :: i

	mag = 0.d0
	Do i=1,3
		mag = mag + (ntruevec(i)**2)
	enddo
	mag = dsqrt(mag)

	ntruevec = ntruevec/mag

!	mag = 0.d0
!	Do i=1,3
!		mag = mag + (ntruevec(i)**2)
!	enddo
!
!	Print*,"mag of nvec=",mag

	RETURN
end subroutine normalise

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!	Subroutine to calculate physical parameters between SPH
!	particle and MCRT packet
!-----------------------------------------------------------------------
subroutine get_coords(xi,ovecdum,l)
	use constants
	use coordinates
	implicit none
		Real*8,dimension(3),intent(in) :: xi
		Real*8,dimension(3),intent(in) :: ovecdum
		Real*8,intent(in) :: l
		integer :: i
		Real*8 :: dum

! Calculate interaction parameters:

!Closest approach, t:
	i=1
	t = 0.d0
	Do i=1,3
		t = t + (xi(i) - ovecdum(i))*ntruevec(i)
	enddo
!Print*,"t=",t
!Impact parameter c:
	i=1
	dum = 0.d0
	Do i=1,3
		dum = dum + (((ovecdum(i)+t*ntruevec(i))-xi(i))**2)
	enddo
	dum = SQRT(dum)/h
	c_dist = min(dum,zeta)
!Print*,"c=",c_dist
!Radial dist to start position:
	i=1
	dum = 0.d0
	Do i=1,3
		dum = dum + ((ovecdum(i)-xi(i))**2)
	enddo
	dum = SQRT(dum)/h
	s0 = min(dum,zeta)
!Print*,"ovec=",ovec(1:3)
!Print*,"xivec=",xivec(1:3)
!Print*,"h=",h
!Print*,"zeta=",zeta
!Print*,"S0=",s0
!STOP

!Radial dist to end position:
	i=1
	dum = 0.d0
	Do i=1,3
		dum = dum + ((ovecdum(i) + l*ntruevec(i)-xi(i))**2)
	enddo
	dum = SQRT(dum)/h
	sl = min(dum,zeta)


!Add sanity check for rounding errors:
	If(sl.lt.c_dist) then
		sl = c_dist
	endif

	If(s0.lt.c_dist) then
		s0 = c_dist
	endif

	RETURN
end subroutine get_coords
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!		Subroutine of total Column Density of ray
!		denoted by sigma
!	Takes contributions to column density from all intersecting particles:
!-----------------------------------------------------------------------
subroutine sig (l,ovecdum,sigma)
	use constants
	use coordinates
	implicit none
		Real*8, intent(in) :: l
		ReaL*8, dimension(3),intent(in) :: ovecdum
		Real*8, intent(out) :: sigma
		Real*8 :: W0, Wl
		Real*8 :: sigmadum
		Logical :: Logicaldum
		Real*8,dimension(3) :: xi
		integer :: i

	sigmadum = 0.d0
	sigma = 0.d0

	Do i=1,nintersects
		xi(1:3) = xi_intersects(1:3,i)
		
	!important that we set new h and m for each xi:
		call set_particle_params(i)
!Print*,"[@sig:] m=",m
!Print*,"[@sig:] h=",h
!Print*,"[@sig:] zeta=",zeta
	
	!Set physical constants:
		call get_coords(xi,ovecdum,l)

!Print*,"[@sig:] nvec=",nvec(1:3)
!Print*,"[@sig:] l=",l
!Print*,"[@sig:] c_dist=",c_dist

!Print*,"[@sig:] s0=",s0
!Print*,"[@sig:] sl=",sl
	!Get column density contributions:	
		call col_dens(s0,W0)
		call col_dens(sl,Wl)
!Print*,"[@sig:] w0=",w0
!Print*,"[@sig:] wl=",wl
	!If contributions on different sides of centre then add,
	! else subtract them to get correct integral:
		If ((t>=0.d0).and.(t<=l)) then
			sigmadum = (m/(h**2)) * (W0 + Wl)
		else
			sigmadum = (m/(h**2)) * (abs(W0 - Wl))
		endif
!Print*,"[@sig:] sigmadum=",sigmadum

		sigma = sigma + sigmadum
!Print*,"[@sig:] sigma=",sigma		
	enddo
!Print*,
!Print*,"[@sig:] sigma=",sigma
!STOP	
!Print*,

	RETURN
end subroutine sig
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!		Subroutine first derivative of Column density
!		As given by O.Lomax & A.P.Whitworth 2016
!-----------------------------------------------------------------------
subroutine sigp (ln,ovecdum,sigmap)
	use constants
	use coordinates
	implicit none
		
		Real*8,intent(in) :: ln					!Last guess of length root
		Real*8,dimension(3),intent(in) :: ovecdum
		Real*8, intent(out) :: sigmap
		Real*8 :: rho
		real*8,dimension(3) :: vec

	vec = (ovecdum + ln*nvec)

	call density (vec,rho)

	sigmap = rho 

!Print*,"[@sigp] sigmap =",sigmap
!Print*,"[@sigp] ln =",ln	

	RETURN
end subroutine sigp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!		Subroutine second derivative of Column density
!		As given by O.Lomax & A.P.Whitworth 2016
!-----------------------------------------------------------------------
subroutine sigpp (ln,lnm1,ovecdum,sigmapp)
	use constants
	use coordinates
	implicit none
		
		Real*8,intent(in) :: ln,lnm1				!Last guess of length root, and impact parameter c
		Real*8,dimension(3),intent(in) :: ovecdum
		Real*8, intent(out) :: sigmapp
		Real*8 :: sigma, sigmam1, sigmap
		Real*8 :: diff

	diff = (lnm1-ln)	

	call sig(ln,ovecdum,sigma)
	call sig(lnm1,ovecdum,sigmam1)
	call sigp(ln,ovecdum,sigmap)

	sigmapp = 2.d0*( ((sigmam1-sigma)/(diff**2)) &
		& - ((sigmap)/diff))


!Print*,"[@sigpp] sigmapp =",sigmapp
!Print*,"[@sigpp] sigma =",sigma
!Print*,"[@sigpp] sigmam1 =",sigmam1
!Print*,"[@sigpp] sigmap =",sigmap
!Print*,"[@sigpp] diff =",diff
!Print*,"[@sigpp] lnm1 =",lnm1
!Print*,"[@sigpp] ln=",ln

	RETURN
end subroutine sigpp

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!		Subroutine of M4 cubix spline kernel derivative
!		As given by O.Lomax & A.P.Whitworth 2016
!-----------------------------------------------------------------------
subroutine kp (s,wp)
	use constants
	use coordinates
	implicit none
		Real*8,intent(in) :: s
		Real*8, intent(out) :: wp
		integer :: i
		Real*8 :: dum


	If (s<=1.d0) then
		wp = ((3.d0/4.d0)*((2.d0-s)**2)) - (3.d0*((1.d0-s)**2))
	else if ((s>1.d0).and.(s<=2.d0)) then
		wp = ((3.d0/4.d0)*((2.d0-s)**2))
	else 
		wp = 0.d0
	endif
	wp = (-1.d0/Pi) * wp

!Print*,"[@kp:] Wp=",wp
	RETURN
end subroutine kp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!		Subroutine of SPH M4 Cubic Spline kernel, denoted by w
!-----------------------------------------------------------------------
subroutine kernel (s, w)
	use constants
	use coordinates
	implicit none
		Real*8, intent(in) :: s
		Real*8, intent(out) :: w

!Print*,"[@kernel] s=",s
	
	If (s<=1.d0) then
		w = ((1.d0/4.d0)*((2.d0-s)**3)) - ((1.d0-s)**3)
	else if ((s>1.d0).and.(s<=2.d0)) then
		w = ((1.d0/4.d0)*((2.d0-s)**3))
	else 
		w = 0.d0
	endif

	w = (1.d0/Pi) * w

!Print*,"[@kernel] w=",w

	RETURN
end subroutine kernel
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!		Subroutine of Column Density, denoted by W
!-----------------------------------------------------------------------
subroutine col_dens (s,W)
	use constants
	use coordinates
	implicit none		
		Real*8,intent(in) :: s					!Last guess of length root in radial distance
		Real*8, intent(out) :: W
		Real*8 :: c2, log1,log2,log3,log4, s2,sqrte		!Pre-evaluated forms for speed of calculation
		Real*8 :: bah
		Real*8, external :: winner, wouter

!Print*,"[@col_dens] log=",log1
!Print*,"[@col_dens] sqrt=",sqrte
!Print*,"[@col_dens] s=",s
!Print*,"[@col_dens] c",c

	If ((c_dist<=1.d0).and.(s<=1.d0)) then
	!Winner(c,s)
		W = winner(c_dist,s)
!Print*,"Inner!"
	else if ((c_dist<=1.d0).and.(s>1.d0)) then
	!Wmixed(c,s)
		W = Winner(c_dist,1.d0) + Wouter(c_dist,s) - Wouter(c_dist,1.d0)
!Print*,"Mixed!"
	else if ((c_dist>1.d0).and.(s>1.d0)) then
	!Wouter(c,s)
		W = wouter(c_dist,s)
!Print*,"Outer!"
	else if (c_dist>2.d0) then
		W = 0.d0
!Print*,"Other!"
	endif

!Print*,"1/pi missing at col_dens!!!!"
	W = (1.d0/Pi)*W

	RETURN
end subroutine col_dens
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!		Function for Inner Column Density
!-----------------------------------------------------------------------
Real*8 function winner(c1,s)
	use constants
	use coordinates
		Real*8,intent(in) :: c1,s
		Real*8 :: c2, log1, s2,sqrte		!Pre-evaluated forms for speed of calculation
		Real*8 :: bah
!Print*,"winner!"
	c2 = c1**2
	s2 = s**2

!Print*,"[@winner:] s=",s
!Print*,"[@winner:] c=",c

	sqrte = DSQRT((s2)-(c2))

	If(isnan(sqrte).eqv..true.)then
		Print*,"[@winner:] sqrte NaN!! Program will abort!"
		Print*,"[@winner:] sqrte=",sqrte
		Print*,"[@winner:] s=",s
		Print*,"[@winner:] c=",c1
		STOP
	endif

!Fix Log issues:
	If (c1==0.d0) then
	!Because log1 terms are multiplied by c, if c==0 this prevents NaN issue!
		log1 = 0.d0
	else
		log1 = Log((s + sqrte)/c1)
	endif

	winner = sqrte*( ((3.d0/16.d0)*(s2*s)) - (0.5d0*s2) &
	& + ((9.d0/32.d0)*(c2*s)) - c2 + 1.d0) &
	& + ((9.d0/32.d0)*(c2**2)*log1)

	RETURN
end function winner
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!		Function for Outer Column Density
!-----------------------------------------------------------------------
Real*8 function wouter(c1,s)
	use constants
	use coordinates
		Real*8,intent(in) :: c1,s
		Real*8 :: c2, log1, s2,sqrte		!Pre-evaluated forms for speed of calculation
		Real*8 :: bah

!Print*,"wouter!"
	c2 = c1**2
	s2 = s**2

	sqrte = DSQRT((s2)-(c2))
	If(isnan(sqrte).eqv..true.)then
		Print*,"[@wouter:] sqrte NaN!! Program will abort!"
		Print*,"[@wouter:] sqrte=",sqrte
		Print*,"[@wouter:] s=",s
		Print*,"[@wouter:] c=",c1
		STOP
	endif

!Fix Log issues:
	If (c1==0.d0) then
	!Because log1 terms are multiplied by c, if c==0 this prevents NaN issue!
		log1 = 0.d0
	else
		log1 = Log((s + sqrte)/c1)
	endif

	wouter =  sqrte*( ((-1.d0)*(1.d0/16.d0)*s2*s) + (0.5d0*s2) &
	& - ((3.d0/32.d0)*(c2 + 16.d0)*s) + c2 + 2.d0) &
	& - ((3.d0/32.d0)*c2*(c2+16.d0)*log1)

	RETURN
end function wouter
