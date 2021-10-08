!Andrew T. Hannington
!ATH4
!120000819
!ath4@st-andrews.ac.uk
!
!Program to use a second order Newton Raphson method on SPH kernel Column Density.
!This allows for the distance a MCRT photon packet should travel for a given optical distance (tau)
!to be calculated.

!-----------------------------------------------------------------------
!
!=======================================================================
!************************ MAIN PROGRAM	********************************
!=======================================================================
subroutine dist_finder (tau,lf,ExitBool)
	use constants
	use coordinates
	implicit none
		Real*8, intent(in) :: tau				!Optical distance cast by MCRT main
		Real*8, intent(out) :: lf				!Distance packet should travel from otruevec
		Logical, intent(out) :: ExitBool			!Has the packet been cast into a vacuum (exited) and thus needs to be killed.
		Real*8 :: ln,lnp1,lnm1					!Current, next, and last guess of distance packet should travel
		Real*8 :: l0						!Initial over-estimate of how far MCRT packet should travel
		integer :: i,j,k,jj					
		Real*8 :: step						!Step iterator should take in distance finder tests
		Real*8 :: x,xp1,xm1					!Sigma for different ln
		Real*8 :: testunit					!Unit being changed in distance finder tests
		Logical :: BracketBool					!Does a solution exist between two distances?
		integer :: endk						!End iterator for distance finder tests
		Real*8 :: l_av						!Average of last 10 guesses of l
		Real*8,dimension(3) :: ovecnew				!New ovec from under estimate subroutine
		Real*8 :: l0new,taudumnew				!New l0 and dummy optical depth from under estimate subroutine
		Real*8 :: tauprime					!Optical depth from origin to a given l
		Real*8 :: ldum						!Dummy l which gives the sum of distances the ovec has been moved in the under-estimate subroutine.
									!This is then added to the distance l found from ovec, to give the total distance from otruevec.

!Print*,"ovec before solution=", otruevec(1:3)

!Initialise:
	l_av = 0.d0
	ldum = 0.d0
	lf = 0.d0

!Set initial taudum to be the optical depth cast in MCRT main.
!This allows for taudum to change but tau to be kept constant for reference
	taudum = tau 
	!Set tau*(chi^(-1))
	tauchim1 = (taudum*(chi**(-1)))

!Initialise exitbool as false
	ExitBool = .FALSE.

!	Open(1,file=filename)

!	Allocate(xi_array(3,n_particles))

!	xi_array(1:3,1) = xi1
!	xi_array(1:3,2) = xi2
!	xi_array(1:3,3) = xi3

!	allocate(xi_intersects(n_particles))
!	xi_intersects = .TRUE.

!	tau = 0.1d0
!	taudum = tau 
!	otruevec = (/0.0d0,1.5d0,0.0d0/)
!	ntruevec = (/0.0d0,-1.d0,0.0d0/)


!	step = 1.d-3
!	zeta = 2.d0*h	


!	endk = 5000

!	do k=1,endk
!		tau = 0.1d0
!		taudum = tau 
!		ldum = 0.d0
		l_av = 0.d0

!		ExitBool = .FALSE.
!		Print*,"Percentage complete=",(dble(k)/dble(endk))*100.d0,"%"

!		testunit = dble(k)*step
!		otruevec(1) = testunit
		!otruevec(1) = testunit

!Print*,"testunit=",testunit



	!Normalise unit vectors:
		call normalise

	!Initialise vectors:
		ovec = otruevec
!		nvec = ntruevec

	!Get initial overestimate of l, l0
		call get_l0 (ovec,l0)


	!Check answer is still physical:
		If (l0.lt.0.d0) then
!			Print*,"l0 < 0! Packet has escaped! Will exit!"
			ldum = 0.d0
			ovec = otruevec
			ExitBool = .TRUE.
			GOTO 22
		else if (l0 .ge. 1.d125) then
!			Print*,"l0 -> +inf ! Packet has escaped! Will exit!"
			ldum = 0.d0
			ovec = otruevec
			ExitBool = .TRUE.
			GOTO 22
		endif


	!Calculate column density from ovec to l0
	! in order to check whether l0 is an underestimate
		call sig(l0,ovec,x)
		tauprime = x * chi

	!Check l0 provides sensible/physical solution
	! i.e. tau' > tau
	! If not, then call under_estimate subroutine
		If (tauprime.lt. tau) then
!			Print*,"Under Est!"
			call under_est (tau,l0,ovec,taudum,l0new,ovecnew,taudumnew,ExitBool,ldum)
			l0 = l0new
			ovec = ovecnew
			taudum = taudumnew
			tauchim1 = (taudum*(chi**(-1)))
		endif

		If(ExitBool .eqv. .TRUE.) then
!			Print*,"under_est found no solution! Packet has escaped! Will exit!"
			GOTO 22
			ovec = otruevec
			ldum = 0.d0
			lf = 0.d0
		endif	

		12 continue

	!Check a solution exists between 0 and l0:
		call root_bracket_test (0.d0,l0,BracketBool)
	!If no solution exists, l0 may be too large. Thus, reduce l0 until root is bracketed
	!else no solution exists if l0 goes too small.
		If (BracketBool .eqv. .FALSE.) then
!		Print*,"Over Est!"
!			Print*,"Solution Lies Outside of Current Bracket!"
!			Print*,"Attempting to force convergence!"
			l0 = l0*convergence_factor

			If(l0.ge.(l0_min))then
				GOTO 12
			else if (l0.lt.(l0_min)) then
!				Print*,"L0 cannot Converge! Packet has escaped! Will Exit!"
				ExitBool = .TRUE.
				lf = 0.d0
				ldum = 0.d0
				ovec = otruevec
				GOTO 22
			endif
		endif		

	!Set tau*(chi^(-1))
	tauchim1 = (taudum*(chi**(-1)))
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
	!Begin Newton-Raphson iteration to find l:
		i=0
		Do 
			i=i+1

		!Set lnm1 and ln based on bracketing of solution:
		!Ensures the solution is always contained between l guesses.
			If (i==1) then
				lnm1 = 0.d0
				ln=l0		
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

		!Some preventative measures to avoid solution entering negative l branch:
			If (ln.lt.(0.d0)) then
				ln = lneg
	!			Print*,
	!			Print*,"***!!!***"
	!			Print*,"Negative l! New ln=",ln
	!			call sig(lnp1,xtest)
	!			Print*,"sigma=",xtest,"tauprime=",tauprime
	!			Print*,"***!!!***"
			endif

		!Check whether solution exists in current range:
			call root_bracket_test (lnm1,ln,BracketBool)
		
			If(BracketBool .eqv. .false.) then
				Print*,"[@Dist_finder:] Warning! NR iteration has diverged from root! Solution cannot be found!"
				Print*,"Program will abort!"
				STOP
			endif


		!!-----------------------------------------------------------------------
		!	Call second order NR calculation subroutine
		!-----------------------------------------------------------------------
			call second_order (ln,lnm1,lnp1)
		!!-----------------------------------------------------------------------


		!More sanity checks:
			If(isnan(lnp1) .eqv. .TRUE.) then
				Print*,"[@Dist_finder:] Lnp1 NaN! Program will abort!"
				STOP
			else if ((lnp1.ge.(1.d125)).or.(lnp1.le.(-1.d125))) then
				Print*,"[@Dist_finder:] Lnp1 inf! Program will abort!"
				STOP
			endif

!			call sig(lnp1,ovec,x)
!			tauprime = x * chi

		!Check for solution convergence:
			If (abs(lnp1-ln)<=tol) then !((abs(tauprime-tau)/tau).le.tol)then!
!				Print*,"Tolerance Reached!"
				GOTO 21
			endif


		!If the number of convergence attempts is too high (a solution cannot reach tolerance of convergence)
		!take the average of the last 10 steps:
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
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

	!This is the solution was found branch:
		21 CONTINUE
	
		lf = lnp1
		ExitBool = .False.

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
	!This is the solution was not found (exitbool === true) branch:
		22 Continue

	!Account for any origin changes from l0 being too small:
		lf = lf + ldum

!Print*,"lnp1=",lnp1
!Print*,"lf=",lf
!Print*,"l0=",l0
!Print*,"ldum=",ldum
!call root_bracket_test (0.d0,(lf*1.01d0),BracketBool)
!Print*,"0 to lf bracketbool=",BracketBool
!call root_bracket_test (0.d0,ldum,BracketBool)
!Print*,"0 to ldum bracketbool=",BracketBool
!Print*,
!Print*,"**!!!**"
!Print*,"ldum=",ldum
!Print*,"lf=",lf
!Print*,"**!!!**"
!Print*,


!		write(1,*) (/testunit,lf/)
		

!	enddo

!	close(1)

!		Print*,"***===***"
!		Print'(A20,ES10.2E3)',"Final l0=",l0
!		Print'(A20,ES10.2E3)',"Final Root Guess=",lf
!		Print*,"ExitBool=",Exitbool
!		Print*,"***===***"
	
!	Print*,"ovec after solution=", ovec(1:3)

	RETURN
end subroutine dist_finder
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!		Subroutine for l1 l2 bracketing root test
! If column densities are of l1 and l2 are on oppiste sides of the root
! the root lies between them and so a solution can be found
! => bracketbool === true. 
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
!		Subroutine for finding new l0 if initial guess is 
!		too small
!-----------------------------------------------------------------------
subroutine under_est (tau,l0old,ovecold,taudumold,l0new,ovecnew,taudumnew,ExitBool,ldum)
	use constants
	use coordinates
	implicit none
		Real*8,intent(in) :: tau
		Real*8, intent(in) :: l0old,taudumold
		Real*8,dimension(3),intent(in) :: ovecold
		Real*8, intent(out) :: l0new,taudumnew
		Real*8,dimension(3),intent(out) :: ovecnew
		Logical,intent(out) :: ExitBool
		Real*8,intent(out) :: ldum
		Integer :: j
		Real*8 :: l0sum
		Real*8 :: tauprime
		Real*8 :: x 
		Logical :: BracketBool,VacuumBool,nointerceptbool
		integer :: i
		Real*8,dimension(3) :: xi
		Real*8 :: tdum,t_error
		Real*8 :: rho		
		Real*8,dimension(3) :: end_position_vector

		taudumnew = taudumold
		l0sum = l0old
		ovecnew = ovecold
		taudumnew = taudumold

		ldum = 0.d0

		j = 0
		11 CONTINUE
		j=j+1

	!Calculate new guess of l0 and sum them (skip over first attempt as l0 already known):
		If(j.ne.(1)) then
			call get_l0 (ovecnew,l0new)
			l0sum = l0sum + l0new
		endif

!Print*,"[@Distance Finder:] tauprime taudum!"
!Print*,"[@Distance Finder:] tauprime taudum attempt=",j

	!Check for infinite loops:
		If(j.ge.(loopmax-10))then

			Print*," "
			Print*,"[@Distance Finder:] tauprime taudum!"
			Print*,"[@Distance Finder:] step j=",j
			Print*,"[@Distance Finder:] l0 =", l0new
			Print*,"[@Distance Finder:] ldum =", ldum
			Print*,"[@Distance Finder:] ovec =", ovecnew(1:3)
			Print*,"[@Distance Finder:] nvec =",nvec(1:3)
			Print*,"[@Distance Finder:] tau =", tau
			Print*,"[@Distance Finder:] taudum =", taudumnew
			Print*,"[@Distance Finder:] tauprime =", tauprime
			Print*,"[@Distance Finder:] xi_intersects=",xi_intersects(1:3,1:nintersects)
			Print*,"[@Distance Finder:] ExitBool=",ExitBool

		call root_bracket_test (0.d0,ldum,BracketBool)

			Print*,"[@Distance Finder:] BracketBool=",BracketBool

		call density (ovecnew,rho)

			Print*,"[@Distance Finder:] Current Location rho=",rho

		call sig(l0new,ovecnew,x)

			Print*,"[@Distance Finder:] Current Location sigma=",x

		call vacuum_check(ovecnew,VacuumBool)

			Print*,"[@Distance Finder:] Current Location VacuumBool=",VacuumBool

			Print*," "
			If(j.ge.loopmax) then
				Print*,"[@Distance Finder:] tauprime taudum!"
				Print*,"Loop Limit exceeded! Infinite Loop Detected! Program will abort!"
				STOP
			endif
		endif

	!Check answer is still physical:
		If (l0new.le.0.d0) then
!			Print*,"Packet has escaped! Exiting Calc!"
			ldum = 0.d0
			ovecnew = otruevec
			ExitBool = .TRUE.
			GOTO 32
		else if (l0new .ge. 1.d125) then
			ldum = 0.d0
			ovecnew = otruevec
			ExitBool = .TRUE.
			GOTO 32
		endif

!Print*,"[@under est] STOP - FIX BELOW HERE!!!"
!STOP
	

	!Setup a vector from ovecnew to end of path at l0new
		end_position_vector  = ovecnew + l0new*ntruevec

	!The above end vector can then be used to check for vacuums
	!and intersects. If in a vacuum, and no intersects remain:
	! there is no solution to be found for this tau and packet
	! should be killed.
		call vacuum_check(end_position_vector,VacuumBool)
!Print*,"VacuumBool=",VacuumBool
		Do i=1,nintersects
			xi(1:3) = xi_intersects(1:3,i)
			call set_particle_params(i)

			call get_t(xi,end_position_vector,tdum)
			call get_coords(xi,end_position_vector,0.d0)

			t_error = DSQRT(((h*zeta)**2)-((h*c_dist)**2))

			If (((tdum + t_error).ge. (0.d0)).and.(c_dist.le.zeta)) then
!Print*,"Intersecting particle i=",i
!Print*,"Intersecting particle xi=",xi(1:3)
!Print*,"Intersecting particle end_position_vector=",end_position_vector(1:3)
!Print*,"Intersecting particle tdum=",tdum
!Print*,"Intersecting particle t_error=",t_error
!Print*,"Intersecting particle (tdum + t_error)=",(tdum + t_error)
!Print*,"Intersecting particle c_dist=",c_dist

				nointerceptbool = .false.
				GOTO 23
			else
				nointerceptbool = .true.
			endif
		enddo

		23 continue
	!Kill if in vacuum and no more intercepts. Packet has escaped!
		If ((VacuumBool .eqv. .true.).and.(nointerceptbool.eqv..true.)) then
!			Print*,"Ovec out of system! l0 infinite! Exiting calc!"
			ldum = 0.d0
			ovecnew = otruevec
			ExitBool = .TRUE.
			GOTO 32
		endif

	!Calculate optical depth from otruevec to total l0:
		call sig(l0sum,otruevec,x)
		tauprime = x * chi

!Print*,"x=",x,"tauchim1=",tauchim1
!Print*,"tau",tau,"tauprime",tauprime

	!Check l0 provides a solution that reaches the cast optical depth
	! i.e. tau' > tau
	! If not, then move origin and re-guess with new origin and taudum
		If (tauprime.lt. tau) then
			ovecnew = ovecnew + l0new*ntruevec
			taudumnew = taudumnew - tauprime
			ldum = ldum + l0new
			GOTO 11
		endif

		32 continue

!		ovec = otruevec
!		l0new = l0sum

	RETURN
end subroutine under_est
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

!Get column density, first and second derivatives
	call sig (ln,ovec,x)
	call sigp (ln,ovec, xp)
	call sigpp (ln,lnm1,ovec,xpp)

!Numerator:
	numer = 2.d0*(x - tauchim1)

!Left and right of denominator SQRT
	left = xp**2
	right = 2.d0*(x - tauchim1)*xpp

!Check SQRT won't be NaN and sigma, sigma', sigma'' aren't infinite:
	If(left.lt.right)then
		Print*," "
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: SQRT NaN! Program will abort!"
		Print*,"***!!!***"
		Print*," "
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*," "
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*," "

		call root_bracket_test (lnm1,ln,BracketBool)
		Print*,"ln=",ln
		Print*,"lnm1=",lnm1
		Print*,"tauchim1=",tauchim1
		Print*,"Root is Bracketed (T/F):",BracketBool 
		STOP
	else if ((x > 1.d125).or.(x < -1.d125)) then
		Print*," "
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: sigma 	INF! Program will abort!"
		Print*,"***!!!***"
		Print*," "
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*," "
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*," "
		STOP

	else if ((xp > 1.d125).or.(xp < -1.d125)) then
		Print*," "
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: sigma prime INF! Program will abort!"
		Print*,"***!!!***"
		Print*," "
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*," "
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*," "
		STOP

	else if((xpp > 1.d125).or.(xpp < -1.d125)) then
		Print*," "
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: sigma prime-prime INF! Program will abort!"
		Print*,"***!!!***"
		Print*," "
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*," "
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*," "
		STOP


	endif	


!Calculate new lnp1:
	sqrte = DSQRT(left - right)

	denom = xp + sqrte

	lnp1 = ln - (numer/denom)


!check new lnp1 is physical:
	If ((lnp1 > 1.d125).or.(lnp1 < -1.d125)) then
		Print*," "
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: lnp1 INF! Program will abort!"
		Print*,"***!!!***"
		Print*," "
		Print*,"lnp1=",lnp1
		Print*,"numerator=",numer
		Print*,"denominator=",denom
		Print*," "
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*," "
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*," "
		STOP
	else if (isnan(lnp1) .eqv. .TRUE.) THEN
		Print*," "
		Print*,"***!!!***"
		Print*,"[@Second Order] WARNING: lnp1 NaN! Program will abort!"
		Print*,"***!!!***"
		Print*," "
		Print*,"lnp1=",lnp1
		Print*,"numerator=",numer
		Print*,"denominator=",denom
		Print*," "
		Print*,"sl=",sl,"ln=",ln,"c=",c_dist
		Print*," "
		Print*,"(x - tauchim1)=",(x - tauchim1)
		Print*,"sigma=",x
		Print*,"sigma p =",xp
		Print*,"sigma pp =",xpp
		Print*,"left =", left
		Print*,"right =", right
		Print*," "
		STOP
	endif

	RETURN
end subroutine second_order
