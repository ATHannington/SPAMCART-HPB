!Andrew T. Hannington
!ATH4
!120000819
!ath4@st-andrews.ac.uk
!
!Created: 20/2/18
!
!Dummy program to use Newton Raphson method on SPH kernel Column Density

!-----------------------------------------------------------------------
!
!=======================================================================
!************************ MAIN PROGRAM	********************************
!=======================================================================
program dist_finder
	use constants
	use coordinates
	implicit none
		Real*8:: tau
		Real*8 :: lf
		Logical :: ExitBool
		Real*8 :: ln,lnp1,lnm1
		Real*8 :: l0
		integer :: i,j,k,jj
		Real*8 :: step
		Real*8 :: x,xp1,xm1
		Real*8 :: testunit
		Logical :: BracketBool
		Logical :: SolBool
		Real*8 :: tdum
		Real*8, dimension(3) ::	xi
		Real*8 :: tauprime							!Effective optical depth of solution
		Real*8 :: ldum,l0sum
		integer :: endk
		Logical :: NoSolBool
		Real*8 :: lprojected
		Real*8 :: rho
		Logical :: vacuumbool
		Real*8 :: taudumsum
		Logical :: nointerceptbool
		Real*8 :: t_error
		Real*8 :: l_av


	ldum = 0.d0
	l_av = 0.d0

!	taudum = tau 

	ExitBool = .FALSE.
	NoSolBool = .FALSE.
	nointerceptbool = .true.

	Open(1,file=filename)

	allocate(xi_array(3,n_particles))
!Setup Particle Parameter arrays:
	allocate(h_array(n_particles))
	allocate(m_array(n_particles))
	allocate(h_intercept_array(n_particles))
	allocate(m_intercept_array(n_particles))

!Setup SPH Particles:	
	xi_array(1:3,1) = xi1(1:3)
!	xi_array(1:3,2) = xi2(1:3)
!	xi_array(1:3,3) = xi3(1:3)
	h_array = hconst
	m_array = mconst

!Setup intersect array:
	allocate(xi_intersects(3,n_particles))

	xi_intersects = xi_array
	nintersects = 1
	h_intercept_array = hconst
	m_intercept_array = mconst


	tau = 0.1d0
	taudum = tau 
	otruevec = (/0.0d0,1.5d0,0.0d0/)
	ntruevec = (/0.0d0,-1.d0,0.0d0/)


	step = 1.d-3
	endk = 25000

	do k=1,endk
		ExitBool = .FALSE.
		NoSolBool = .FALSE.
		nointerceptbool = .true.

		tau = 0.1d0
		taudum = tau 
		ldum = 0.d0
		l_av = 0.d0

		Print*,"Percentage complete=",(dble(k)/dble(endk))*100.d0,"%"

		testunit = dble(k)*step
!		mconst = testunit
!		m_array = mconst
!		m_intercept_array = mconst
		taudum = testunit
!		otruevec(1) = testunit
		!otruevec(1) = testunit

!Print*,"testunit=",testunit
		j=0

		!Normalise unit vectors:
		call normalise

		!Initialise vectors:
		ovec = otruevec
		nvec = ntruevec

		tauchim1 = (taudum*(chi**(-1)))

		taudumsum = 0.d0
		ldum = 0.d0
		lf = 0.d0
		l0 = 0.d0
		l0sum = 0.d0

		11 CONTINUE
		j=j+1
!Print*,"[@Distance Finder:] tauprime taudum!"
!Print*,"[@Distance Finder:] tauprime taudum attempt=",j
		If(j.ge.(loopmax-10))then

			Print*,
			Print*,"[@Distance Finder:] tauprime taudum!"
			Print*,"[@Distance Finder:] step j=",j
			Print*,"[@Distance Finder:] l0 =", l0
			Print*,"[@Distance Finder:] ldum =", ldum
			Print*,"[@Distance Finder:] ovec =", ovec(1:3)
			Print*,"[@Distance Finder:] nvec =",nvec(1:3)
			Print*,"[@Distance Finder:] tau =", tau
			Print*,"[@Distance Finder:] taudum =", taudum
			Print*,"[@Distance Finder:] taudumsum =", taudumsum
			Print*,"[@Distance Finder:] tauprime =", tauprime
			Print*,"[@Distance Finder:] xi_intersects=",xi_intersects(1:3,1:n_particles)
			Print*,"[@Distance Finder:] ExitBool=",ExitBool

		call root_bracket_test (0.d0,ldum,BracketBool)

			Print*,"[@Distance Finder:] BracketBool=",BracketBool

		call density (ovec,rho)

			Print*,"[@Distance Finder:] Current Location rho=",rho

		call sig(l0,ovec,x)

			Print*,"[@Distance Finder:] Current Location sigma=",x

		call vacuum_check(ovec,VacuumBool)

			Print*,"[@Distance Finder:] Current Location VacuumBool=",VacuumBool

			Print*,
			If(j.ge.loopmax) then
				Print*,"[@Distance Finder:] tauprime taudum!"
				Print*,"Loop Limit exceeded! Infinite Loop Detected! Program will abort!"
				STOP
			endif
		endif


		!Calculate first guess at distance travelled:
		call get_l0 (ovec,l0)
		l0sum = l0sum + l0
!Print*,"l0=",l0

		12 CONTINUE

	!Print*,
	!Print*,"***+++***"
	!Print*,"[@ dist_finder:]"
	!Print*,"xivec=",xivec(1:3)
	!Print*,"ovec=",ovec(1:3)
	!Print*,"nvec=",nvec(1:3)
	!Print*,"s0=",s0
	!Print*,"c=",c_dist
	!Print*,"t=",t
	!Print*,"l0=",l0
	!Print*,"***+++***"
	!Print*,

		!Check answer is still physical:
		If (l0.lt.0.d0) then
!			Print*,"Packet has escaped! Exiting Calc!"
			lnp1 = -2.d-1
			ldum = 0.d0
			ovec = otruevec
			ExitBool = .TRUE.
			GOTO 22
		else if (l0 .ge. 1.d125) then
!			Print*,"l0 infinite! Exiting calc!"
			lf = 0.d0
			ldum = 0.d0
			ovec = otruevec
			ExitBool = .TRUE.
			GOTO 22
		endif

		call vacuum_check(ovec,VacuumBool)

		Do i=1,nintersects
			xi(1:3) = xi_intersects(1:3,i)
			call set_particle_params(i)

			call get_t(xi,ovec,tdum)
			call get_coords(xi,ovec,0.d0)

			t_error = DSQRT((zeta**2)-(c_dist**2))

			tdum = tdum + t_error				

			If (tdum .ge. (0.d0)) then
				nointerceptbool = .false.
				GOTO 23
			else
				nointerceptbool = .true.
			endif
		enddo

		23 continue

		If ((VacuumBool .eqv. .true.).and.(nointerceptbool.eqv..true.)) then
!			Print*,"Ovec out of system! l0 infinite! Exiting calc!"
			lf = 0.d0
			ldum = 0.d0
			ovec = otruevec
			ExitBool = .TRUE.
			GOTO 22
		endif

		If (UnderEstBool .eqv. .TRUE.) then

			call sig(l0sum,ovec,x)
			tauprime = x * chi

			taudumsum = taudumsum + tauprime

	!Print*,"x=",x,"tauchim1=",tauchim1
	!Print*,"tau",tau,"tauprime",tauprime

			!Check l0 provides sensible/physical solution
			! i.e. tau' > tau
			! If not, then move origin and re-guess with new origin and tau
			If (taudumsum.lt. tau) then
				ovec = ovec + l0*nvec
				taudum = taudum - tauprime

!Print*,
!Print*,"***"
!Print*,"tauprime taudum!"
!Print*,"chi=",chi
!Print*,"tauprime=",tauprime
!Print*,"l0=",l0
!Print*,"ldum=",ldum
!Print*,"***"
!Print*,
				ldum = ldum + l0
				GOTO 11
			endif
		endif


		tauchim1 = (taudum*(chi**(-1)))

		i=0
		Do 
			i=i+1
!Print*,"[@Distance Finder:] convergence!"
!Print*,"[@Distance Finder:] convergence attempt=",i
	!		Print*,
	!		Print*,"Current Iteration=",i

			If (i==1) then
				lnm1 = 0.d0
				ln=l0sum		
			else
				If (lnm1SideBool .eqv. .TRUE.)then
					call root_bracket_test (ln,lnp1,BracketBool)
					If (BracketBool .eqv. .TRUE.) then
						lnm1 = ln
					endif
				else
					lnm1 = ln
				endif
				ln = lnp1
			endif

		!Some critical case preventative measures:
			If (ln.lt.(0.d0)) then
				ln = lneg
	!			Print*,
	!			Print*,"***!!!***"
	!			Print*,"Negative l! New ln=",ln
	!			call sig(lnp1,xtest)
	!			Print*,"sigma=",xtest,"tauprime=",tauprime
	!			Print*,"***!!!***"
			endif

		!Check if packet has no physical positive distance solution
		!and thus has escaped system:
			call root_bracket_test (lnm1,ln,BracketBool)
			If (BracketBool .eqv. .FALSE.) then
	!			Print*,"Solution Lies Outside of Current Bracket!"
	!			Print*,"Attempting to force convergence!"
				l0 = l0*convergence_factor

				If(l0.ge.(l0_min))then
					GOTO 12
				else if (l0.lt.(l0_min)) then
	!				Print*,"L0 cannot Converge! Packet has escaped! Will Exit!"
					ExitBool = .TRUE.
					lnp1 = 0.d0
					ldum = 0.d0
					ovec = otruevec
					GOTO 22
				endif
			endif

			call second_order (ln,lnm1,lnp1)

			If(isnan(lnp1) .eqv. .TRUE.) then
				Print*,"[@Dist_finder:] GUESS IS UNPHYSICAL!!!"
				STOP
			endif


			If (abs(lnp1-ln)<=tol) then
	!			Print*,"Tolerance Reached!"
				GOTO 21
			endif



			If(i.ge.(convergeattmeptlim-10))then
				l_av = l_av + lnp1
!				Print*,
!				Print*,"[@Distance Finder:] step i=",i
!				Print*,"[@Distance Finder:] lnp1 =", lnp1
!				Print*,"[@Distance Finder:] ovec =", ovec(1:3)
!				Print*,"[@Distance Finder:] nvec =",nvec(1:3)
!				Print*,"[@Distance Finder:] taudum =", taudum
!				Print*,"[@Distance Finder:] xi_intersects=",xi_intersects(1:n_particles)
!				Print*,"[@Distance Finder:] ExitBool=",ExitBool
!				Print*,
				If(i.ge.convergeattmeptlim) then
					lnp1 = l_av/10.d0
!					Print*,"[@Distance Finder:]"
!					Print*,"Convergence Limit exceeded!"
				goto 21
				endif
			endif



		enddo

		21 CONTINUE

		lf = lnp1
		ExitBool = .False.
		ExcludeBool = .FALSE.

!		Print*,
!		Print*,"***"
!
!		Print*,"***"
!		Print*,"Final Iterations=",i
!		Print*,"l0=",l0
!		Print*,"Attempts at convergence: j=", j
!		Print*,"***"
!		Print*,
!		Print*,"TestUnit =", TestUnit
!		Print*,
!		Print*,"***===***"
!		Print'(A20,ES10.2E3)',"Final Root Guess=",lf
!		Print*,"***===***"
!		Print*,"***"
!
		22 Continue

	!Account for any origin changes from l0 being too small:
		lf = lf + ldum


!Print*,
!Print*,"**!!!**"
!Print*,"ldum=",ldum
!Print*,"lf=",lf
!Print*,"**!!!**"
!Print*,


		write(1,*) (/testunit,lf/)
		

	enddo

	close(1)

!		Print*,"***===***"
!		Print'(A20,ES10.2E3)',"Final Root Guess=",lf
!		Print*,"ExitBool=",Exitbool
!		Print*,"***===***"

	RETURN
end program dist_finder
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!		Subroutine for l1 l2 bracketing root test
!-----------------------------------------------------------------------
subroutine root_bracket_test (l1,l2,BracketBool)
	use constants
	use coordinates
	implicit none
		Real*8, intent(in) :: l1,l2
		Logical, intent(out) :: BracketBool
		Real*8 :: x1, x2

	call sig(l1,ovec,x1)
	call sig(l2,ovec,x2)

	x1 = (x1 - tauchim1)
	x2 = (x2 - tauchim1)

	If ( ((x1<=0.d0).and.(x2>=0.d0)).or.((x1>=0.d0).and.(x2<=0.d0)) ) then
		BracketBool = .TRUE.
!Print*,"Bracket TRUE"
!Print*,"l1=",l1
!Print*,"l2=",l2
!Print*,"x1=",x1
!Print*,"x2=",x2
	else if ( ((x1<=0.d0).and.(x2<=0.d0)).or.((x1>=0.d0).and.(x2>=0.d0)) ) then
		BracketBool = .FALSE.
!Print*,"Bracket FALSE"
!Print*,"l1=",l1
!Print*,"l2=",l2
!Print*,"x1=",x1
!Print*,"x2=",x2
	endif

	RETURN
end subroutine root_bracket_test
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!		Subroutine for Second Order Newton Raphson Solver
!-----------------------------------------------------------------------
subroutine second_order (ln,lnm1,lnp1)
	use constants
	use coordinates
	implicit none
		Real*8, intent(in) :: ln,lnm1
		Real*8, intent(out) :: lnp1
		Real*8 :: x,xp,xpp
		Real*8 :: denom
		Real*8 :: numer
		Real*8 :: left,right,sqrte
		Logical :: BracketBool

	call sig (ln,ovec,x)
	call sigp (ln,ovec, xp)
	call sigpp (ln,lnm1,ovec,xpp)

	numer = 2.d0*(x - tauchim1)
	
	left = xp**2
	right = 2.d0*(x - tauchim1)*xpp

	If(left.lt.right)then
		Print*,
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: SQRT NaN! Program will abort!"
		Print*,"***!!!***"
		Print*,
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*,
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*,

call root_bracket_test (lnm1,ln,BracketBool)
		Print*,"ln=",ln
		Print*,"lnm1=",lnm1
		Print*,"tauchim1=",tauchim1
		Print*,"Root is Bracketed (T/F):",BracketBool 
		STOP
	else if ((x > 1.d125).or.(x < -1.d125)) then
		Print*,
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: sigma 	INF! Program will abort!"
		Print*,"***!!!***"
		Print*,
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*,
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*,
		STOP

	else if ((xp > 1.d125).or.(xp < -1.d125)) then
		Print*,
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: sigma prime INF! Program will abort!"
		Print*,"***!!!***"
		Print*,
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*,
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*,
		STOP

	else if((xpp > 1.d125).or.(xpp < -1.d125)) then
		Print*,
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: sigma prime-prime INF! Program will abort!"
		Print*,"***!!!***"
		Print*,
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*,
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*,
		STOP


	endif	

	sqrte = DSQRT(left - right)

	denom = xp + sqrte

	lnp1 = ln - (numer/denom)

	If ((lnp1 > 1.d125).or.(lnp1 < -1.d125)) then
		Print*,
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: lnp1 INF! Program will abort!"
		Print*,"***!!!***"
		Print*,
		Print*,"lnp1=",lnp1
		Print*,"numerator=",numer
		Print*,"denominator=",denom
		Print*,
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*,
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*,
		STOP
	else if (isnan(lnp1) .eqv. .TRUE.) THEN
		Print*,
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: lnp1 NaN! Program will abort!"
		Print*,"***!!!***"
		Print*,
		Print*,"lnp1=",lnp1
		Print*,"numerator=",numer
		Print*,"denominator=",denom
		Print*,
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*,
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*,
		STOP
	endif
	RETURN
end subroutine second_order
