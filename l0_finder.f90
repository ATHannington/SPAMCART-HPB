!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!		Subroutine of inverse mean free path 
!-----------------------------------------------------------------------
subroutine inverse_mfp(vec,mfp_m1)
	use constants
	use coordinates
		Real*8,dimension(3),intent(in) :: vec
		Real*8 :: rho
		Real*8, intent(out) :: mfp_m1

	Call density(vec,rho)

	mfp_m1 = chi*rho

!Print*,"[@inverse mfp:] vec=",vec(1:3)
!Print*,"[@inverse mfp:] rho=",rho
!Print*,"[@inverse mfp:] chi=",chi
!Print*,"[@inverse mfp:] mfp_m1=",mfp_m1

	RETURN
end subroutine  inverse_mfp
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!		Subroutine of gradient of inverse mean free path
!-----------------------------------------------------------------------
subroutine grad_inverse_mfp(vec,g_mfp_m1)
	use constants
	use coordinates
		Real*8,dimension(3),intent(in) :: vec
		Real*8,dimension(3), intent(out) :: g_mfp_m1
		Real*8 :: rho, wpdum

		real*8,dimension(3) :: xi,svec
		integer :: i,j
		Real*8 :: dum
		Real*8 :: s
		Logical :: Logicaldum
		Real*8,dimension(3) :: outdum

	outdum = 0.d0
	g_mfp_m1 = 0.d0

!Print*,"[@g_mfp:] vec=",vec(1:3)

	Do i=1,nintersects
		xi(1:3) = xi_intersects(1:3,i)
	!important that we set new h and m for each xi:
		call set_particle_params(i)

		svec = (vec - xi)
!Print*,"[@g_mfp:] svec=",svec(1:3)

		dum = 0.d0
		Do j=1,3
			dum = dum + svec(j)**2
		enddo
		s = DSQRT(dum)/h
!Print*,"[@g_mfp:] s=",s
!		!sets all coordinates apart from sl which is not needed:
		call get_coords(xi,vec,0.d0)

		call kp(s,wpdum)
!Print*,"[@g_mfp:] wpdum=",wpdum

	!s>=0 (radial dist). Avoid NaNs when s=0:
		If (s.le.(0.d0)) then
			outdum = 0.d0
		else
			outdum = chi*(m/h**3)*(1.d0/s)*svec*wpdum
		endif

		g_mfp_m1 = g_mfp_m1 + outdum
		
	enddo

!Print*,"[@g_mfp:] g_mfp_m1=",g_mfp_m1(1:3)

	RETURN
end subroutine  grad_inverse_mfp

!-----------------------------------------------------------------------
!		Subroutine to get l0 as per Lomax eqn (20)
!-----------------------------------------------------------------------
subroutine get_l0 (ovecdum,l0)
	use constants
	use coordinates
	implicit none
		Real*8,dimension(3),intent(in) :: ovecdum
		Real*8,intent(out) :: l0
		Real*8,dimension(3) :: g_mfp_m1
		Real*8 :: mfp_m1, dot
		Real*8 :: delta
		Integer :: i
		Real*8 :: dum

	call inverse_mfp(ovecdum,mfp_m1)
	call grad_inverse_mfp(ovecdum,g_mfp_m1)

!Dot product between ntruevec and gradient of mean free path:
	dot=0.d0
	Do i=1,3
		dot = dot + g_mfp_m1(i)*ntruevec(i)
	enddo

	delta = max( ((mfp_m1**2) + 2.d0*taudum*dot),0.d0)
	
	l0 = (2.d0*overest*taudum)/(mfp_m1 + DSQRT(delta))

!Print*,
!Print*,"[@get_l0] (mfp_m1 + DSQRT(delta))",(mfp_m1 + DSQRT(delta))
!Print*,"[@get_l0] (2.d0*overest*taudum)",(2.d0*overest*taudum)
!Print*,"[@get_l0] ovec=", ovec(1:3)
!Print*,"[@get_l0] nvec",nvec(1:3)
!Print*,"[@get_l0] tau",taudum
!Print*,"[@get_l0] l0",l0
!Print*,"[@get_l0] mfp_m1",mfp_m1
!Print*,"[@get_l0] g_mfp_m1",g_mfp_m1
!Print*,"[@get_l0] dot",dot
!Print*,"[@get_l0] delta",delta
	RETURN
end subroutine get_l0


