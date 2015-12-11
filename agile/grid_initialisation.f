!=======================================================================
!
!     AGILE: initial grid point allocation
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

        subroutine grid_initialisation(old,new,status)

        use agile_parameter_module
        use input_file_module
        use state_vector_module
        use adaptive_grid_module
        use newton_module

        implicit none

        integer, intent(out) :: status
        type(state_vector), intent(inout) :: old
        type(state_vector), intent(out) :: new

!-----------------------------------------------------------------------
!
!     input:
!     old   state vector with initial profiles, no staggering assumed.
!
!     output:
!     adaptive_linint_template   a copy of old state for interpolations
!     new   same state vector consistent with adaptive grid equation
!     status   is set to 1 if the points could not be allocated
!     old   last converged state if status=1
!
!-----------------------------------------------------------------------

        logical :: smooth
        integer :: nq,ny,staggering,iq,ix,iql,iqm,iqu,itmax,it
        real :: dq,q0,f,slowc,conv,dt,tau
        real, dimension(2) :: bpsi
        real, dimension(nnq) :: q,y,psi

!.....read data from input file.........................................
        if (input_action%initial) call input_read_action
        if (input_state%initial) call input_read_state
        nq = input_state%nq
        if (nq.gt.old%nq) then
          if (input_action%message.ge.1)
     1      write(6,*) 'Error: nq too large in grid_initialisation()'
          stop
        endif
        ny = input_state%ny
        if (ny.ne.old%ny) then
          if (input_action%message.ge.1) write(6,11)
11        format('Error: wrong # of variables in grid_initialisation!')
          stop
        endif
        if (any(input_state%staggering(1:ny).eq.-1)) then
          staggering = -1
        elseif (any(input_state%staggering(1:ny).eq.1)) then
          staggering = 1
        else
          staggering = 0
        endif
        if (input_allocation%initial) call input_read_allocation
        tau = input_allocation%tau

!-----create a guess without grid smoothing-----------------------------

!.....copy some parameters into new.....................................
        new%nq = nq
        new%ny = ny
        new%jy(1:nq,1:ny) = old%jy(1:nq,1:ny)
        new%t = old%t

!.....calculate q at points x by a summation of delta(q)................
        smooth = .false.
        call state_scaling(old)
        if ( any(old%yscal(1:old%nq*old%ny).le.0.)
     1    .or. any(old%yscl(1:old%ny).le.0.) ) then
      write(6,*) 'Error: old state not scaled in grid_initialisation()!'
          stop
        endif
        call adaptive_deltaq(smooth,old,old,psi,bpsi,status)
        if (status.ne.0) return
        q(1) = 0.
        do iq=2,old%nq
          q(iq) = q(iq-1) + 1./psi(iq-1)
        enddo

!.....set equal distances in q..........................................
	if (staggering.eq.-1) then
          dq = q(old%nq)/(float(nq)-1.5)
          q0 = 0.
        elseif (staggering.eq.1) then
          dq = q(old%nq)/(float(nq)-1.5)
          q0 = -0.5*dq
	else
          dq = q(old%nq)/float(nq-1)
          q0 = 0.
	endif
        call state_sumdx(old,y)
	new%y(new%jy(1,1)) = y(1)
	do ix=2,nq-1
          iql = 1
          iqu = old%nq
          do while (iqu-iql.gt.1)
            iqm = (iqu+iql)/2
            if ((dq.gt.0.).eqv.(q0+(ix-1)*dq.gt.q(iqm))) then
              iql = iqm
            else
              iqu = iqm
            endif
          enddo
          f = (q0+(ix-1)*dq-q(iql))/(q(iql+1)-q(iql))
	  new%y(new%jy(ix,1)) = y(iql) + f*(y(iql+1)-y(iql))
        enddo
        new%y(new%jy(nq,1)) = y(old%nq)
        new%dx = 0.

!.....interpolate y's at x..............................................
        call adaptive_linint(new,status)
        if (status.ne.0) return

!.....update old state (remains saved in adaptive_linint_template)......
        old = new

!-----relax with grid smoothing-----------------------------------------

!.....initialise relax..................................................
	itmax = 20
	slowc = 0.1d-0
	conv = 0.1d-8
        newton_flag%block = .true.
        newton_flag%analyze = .true.
        newton_flag%banddiag = .true.
        newton_flag%bigrad = .false.
        dt = tau

!.....if (tau.ne.0) do iterations using time-smoothing..................
	do it=1,100
	  if (input_action%message.ge.1)
     1      write(6,*) 'Info: grid point allocation, iteration ',it

!.....relax.............................................................
          status = 0
          new%t = old%t + dt
          call state_scaling(old)
	  call newton_raphson(itmax,slowc,conv,old,
     1      new,status,adaptive_relax)

!.....change timestep...................................................
	  if (status.ne.0) then
	    if (input_action%message.ge.2)
     1        write(6,*) 'Info: no convergence in newton_raphson'
            new = old
            dt = 0.5*dt
	  else
            old = new
            dt = 1.5*dt
	  endif

!.....finished..........................................................
	  if ((dt.eq.0.).or.(dt.lt.tau/1000.).or.(dt.gt.1000.*tau)) then
            if (input_action%message.ge.1) then
	      if (status.ne.0) then
	        write(6,*) 'Warning: relax on allocation guess failed'
	        if (it.eq.1) write(6,*)
     1            'try again with (tau.gt.0.) in section ALLOCATE...'
	        write(6,*) 'Info: writing unrelaxed guess instead'
	      else
	        write(6,*) 'Info: successful Allocation'
              endif
	    endif
            return
	  endif
	enddo

        if (input_action%message.ge.1) then
	  write(6,*) 'Warning: itmax exceeded in alloc.f'
	  write(6,*) 'Info: writing unrelaxed guess instead'
        endif
	status = 1

	end subroutine grid_initialisation

!=======================================================================
