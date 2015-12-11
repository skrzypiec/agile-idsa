!***********************************************************************
!
!     AGILE: time_step_module
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module time_step_module

      use agile_parameter_module
      use input_file_module
      use state_vector_module
      use newton_module

      implicit none

      real, parameter :: time_red=0.95 !time step reduction factor
      real, parameter :: time_inc=1.05 !time step increase factor

      contains

!=======================================================================
!
!     AGILE: time_implicit_step
!
!=======================================================================

	subroutine time_implicit_step(hido,old,new,
     1    hnext,equations_step)

        use eos_module, only: eionglobal

	external equations_step

        real, intent(inout) :: hido
        real, intent(out) :: hnext
        type(state_vector), intent(inout) :: old
        type(state_vector), intent(out) :: new

!-----------------------------------------------------------------------
!
!     input:
!     hido    preferred time step
!     old     state vector with initial state
!     equations_step   routine with equations describing the time step
!
!     output:
!     hido   successful time step
!     new    state vector after hido
!     hnext   proposed next time step
!
!     Note: old has intent(out) because state_write might set
!     numbers that are close to zero to zero exactly.
!
!-----------------------------------------------------------------------

        logical, save :: initial=.true.
        integer :: repeat
        integer, save :: iq=0,iy=0
        integer :: n,imax,nq
        integer, dimension(1) :: iloc
        real, dimension(nnc) :: err

!.....initialise........................................................
        n = old%nq*old%ny
        if (initial) then
          if (input_step%initial) call input_read_step
          if (input_jacobi%initial) call input_read_jacobi
          newton_flag%block = .true.
          newton_flag%analyze = .true.
          newton_flag%banddiag = .true.
          newton_flag%bigrad = .false.
          initial = .false.
        endif

!.....try steps until success...........................................
        repeat = 1
        do while (repeat.ne.0)

!.....display time step.................................................
          if (input_action%message.ge.1) write(6,22) hido,iq,iy

!.....make a guess......................................................
          new = old
          new%t = new%t + hido

!.....take a step.......................................................
          repeat = 0
          call newton_raphson(input_jacobi%itmax,input_jacobi%slowc,
     1      input_jacobi%conv,old,new,repeat,equations_step)

!.....if this time step caused a problem --> reduce it..................
          if (repeat.ne.0) then
            if (input_action%message.ge.1) write(6,33)
            newton_flag%block = .true.
            newton_flag%analyze = .true.
            hido = hido/2.
            if (hido.lt.input_step%minimum) then
              write(6,11)
              call state_write(-1,old)
              stop
            endif

!.....pause.............................................................
            if (input_action%message.ge.4) then
              write(6,*) '(press return to continue)'
              read(5,*)
            endif

!.....calculate the change in the variables.............................
          else
            nq = old%nq
            err(1:n) = abs(new%y(1:n) - old%y(1:n))/old%yscal(1:n)
            err(old%jy(1:nq,1)) = 0. !ignore grid displacements
            
!.....consider rather energy change than temperature change.............
            err(old%jy(1:nq,6)) = abs(new%bout%e(1:nq)-old%bout%e(1:nq))
     1        / (old%bout%e(1:nq) - eionglobal/units%c**2)
            
!.....detect maximum variable change....................................
            iloc = maxloc(err(1:n))
            imax = iloc(1)
            iy = mod(imax,old%ny)
            iq = 1+ (imax-iy)/old%ny
      
!.....adjust time step to variable changes..............................
            if (err(imax).gt.input_step%eps) then
              hnext = max(time_red*hido,input_step%minimum)
            elseif (err(imax).lt.0.5*input_step%eps) then
              hnext = min(time_inc*hido,input_step%maximum)
            else
              hnext = hido
            endif
            write(6,44) 'AGILE changes: ',err(imax)
          endif
        enddo
        initial = .false.

!.....format statements.................................................
11      format(1x,'Error: step size underflow in time_implicit_step!')
22      format(1x,'Info: h_try = ',g10.3,' y(',i3,',',i3,')')
33      format(1x,'Warning: failed time step')
44      format(a24,g12.4)
66      format(1x,'Warning: inaccurate time step in time_step_module')

	end subroutine time_implicit_step

!=======================================================================

      end module time_step_module

!***********************************************************************
