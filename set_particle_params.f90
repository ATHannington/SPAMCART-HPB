!-----------------------------------------------------------------------------------
!Subroutine to setup physical constants for particle being used of particles in intercepts arrays:
!-----------------------------------------------------------------------------------
subroutine set_particle_params(i)
	use constants
	use coordinates
	implicit none
		integer, intent(in) :: i

!	m = m_array(i)
!	h = h_array(i)
	h = h_intercept_array(i)
	m = m_intercept_array(i)
!Print*,"[@particle_params:] mass=",m	

!Sanity checks:
	If (m.le.(0.d0))then
		Print*,"[@Set_particle_params:] Zero mass detected! Program will abort!"
		STOP
	endif

	If (h.le.(0.d0))then
		Print*,"[@Set_particle_params:] Zero smoothing length (h) detected! Program will abort!"
		STOP
	endif

	RETURN
end subroutine set_particle_params
