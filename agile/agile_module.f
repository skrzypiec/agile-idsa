!***********************************************************************
!
!     AGILE: agile_module, provides module-access to agile
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!     AGILE (originally named Adaptive Grid Implicit Leap Extrapolation)
!     has been developed by M. Liebendoerfer at the University of Basel,
!     funded by the Swiss National Science Foundation. It is documented
!     in the Astrophysical Journal 141 (2002) p. 229.
!
!     Second order TVD advection was introduced at CITA, University of
!     Toronto (2003). A new handling of the independent variable
!     'enclosed mass' permits a much improved numerical stability of
!     zones with very small mass content. This improvement is summarised
!     in Fischer et al., Astronomy & Astrophysics 517 (2010) p. 80.
!
!   Agile-IDSA is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License version 3
!   as published by the Free Software Foundation.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!***********************************************************************

      module agile_module

      use agile_parameter_module
      use input_file_module
      use species_module
      use state_vector_module

      implicit none

      external equations_relax,equations_step

!=======================================================================
!
!     AGILE: static memory
!
!=======================================================================

!.....state vector......................................................
      type(state_vector) :: state

!.....time step management..............................................
      real, parameter :: agile_red=0.95 !time step reduction factor
      real, parameter :: agile_inc=1.05 !time step increase factor

!.....exchange terms between AGILE and IDSA.............................
      integer :: istep=0 !current step number
      real, dimension(nnq) :: vstress=0.,venergy=0.,vye=0.
      real, dimension(nut,nnq) :: vynu=0.,vznu=0.
      real, dimension(nnq) :: jrad=0.,hrad=0.,krad=0.

      contains

!=======================================================================
!
!     AGILE: write state vector into file with number "cycle"
!
!=======================================================================

        subroutine agile_write(cycle)

        integer, intent(in) :: cycle

!-----------------------------------------------------------------------

        call state_write(cycle,state)

        end subroutine agile_write

!=======================================================================
!
!     AGILE: initialisation
!
!=======================================================================

        subroutine agile_initialisation

!-----------------------------------------------------------------------

!.....initialise AGILE..................................................
        call agile_start(state,equations_relax)

        end subroutine agile_initialisation

!=======================================================================
!
!     AGILE: time step
!
!=======================================================================

        subroutine agile_step(dt,old,new,dtnext,equations_step)

        use eos_module, only: eionglobal
        use newton_module

        external equations_step

        real, intent(inout) :: dt
        real, intent(out) :: dtnext
        type(state_vector), intent(inout) :: old
        type(state_vector), intent(out) :: new

!-----------------------------------------------------------------------
!
!     input:
!     dt    preferred time step
!     old     state vector with initial state
!     equations_step   routine with equations describing the time step
!
!     output:
!     dt   successful time step
!     new    state vector after dt
!     dtnext   proposed next time step
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
          if (input_action%message.ge.1) write(6,22) dt,iq,iy

!.....make a guess......................................................
          new = old
          new%t = new%t + dt

!.....take a step.......................................................
          repeat = 0
          call newton_raphson(input_jacobi%itmax,input_jacobi%slowc,
     1      input_jacobi%conv,old,new,repeat,equations_step)

!.....if this time step caused a problem --> reduce it..................
          if (repeat.ne.0) then
            if (input_action%message.ge.1) write(6,33)
            newton_flag%block = .true.
            newton_flag%analyze = .true.
            dt = dt/2.
            if (dt.lt.input_step%minimum) then
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
              dtnext = max(agile_red*dt,input_step%minimum)
            elseif (err(imax).lt.0.5*input_step%eps) then
              dtnext = min(agile_inc*dt,input_step%maximum)
            else
              dtnext = dt
            endif
            write(6,44) 'AGILE changes: ',err(imax)
          endif
        enddo
        initial = .false.

!.....format statements.................................................
11      format(1x,'Error: step size underflow in agile_step!')
22      format(1x,'Info: dt_try = ',g10.3,' y(',i3,',',i3,')')
33      format(1x,'Warning: failed time step')
44      format(a24,g12.4)
66      format(1x,'Warning: inaccurate time step in agile_step')

        end subroutine agile_step

!=======================================================================
!
!     AGILE-IDSA interface
!
!=======================================================================

        subroutine agile_idsa(dt,dtnext)

        use delept_module
        use egroup_module
        use spectrum_module, only: src
        use units_module

        real, intent(in) :: dt
        real, intent(out) :: dtnext

!-----------------------------------------------------------------------
!
!       Input:
!       dt ... neutrino time step
!       state ... state vector with up-to-date state_fill_zones
!
!       Output:
!       jrad(i) = 3.*pnu(i)/d(i)                             ![cm^2/s^2]
!       hrad(i) = LE(i)/(4.*units%pi*ar(i)**2*d(i))          ![cm^3/s^3]
!       krad(i) = pnu(i)                                     ![g/cm/s^2]
!       vye(i) = ydot(i)*state%bda(iq)                            ![g/s]
!       venergy(i) = edot(i)*state%bda(iq)                      ![erg/s]
!       vstress(iq) = vdot(i)*state%ada(iq)                  ![g*cm/s^2]
!       dtnext ... proposed next time step                          ![s]
!
!-----------------------------------------------------------------------

        integer :: n,iter,status,i,iq,it,ie,iray
        real :: tmp,zold,znew
        real, dimension(ne) :: ps,rate
        real, dimension(:), allocatable :: ar,dV,d,s,ye,pnu,vdot
        real, dimension(:), allocatable :: sout,edot,yout,ydot
        real, dimension(:,:), allocatable :: ynu,znu,ynudot,znudot
        real, dimension(nut) :: dLmax
        real, dimension(ne,nut) :: rnu
        real :: chg

!.....initialisation....................................................
        if (egroup_empty) call egroup
        n = state%nq-1

!.....convert agile quantities to idsa quantities.......................
        allocate(ar(n),dV(n),d(n),s(n),ye(n),ynu(nut,n),znu(nut,n))
        do i=1,n
          iq = i+1
          ar(i) = state%ar(iq)                                     ![cm]
          if (i.eq.1) then
            dV(i) = 4./3.*units%pi*ar(i)**3                      ![cm^3]
          else
            dV(i) = (ar(i) - ar(i-1))                            ![cm^3]
     &        * 4./3.*units%pi*(ar(i)**2+ar(i)*ar(i-1)+ar(i-1)**2)
          endif
          d(i) = state%bin%n(iq)                               ![g/cm^3]
          s(i) = state%bout%s(iq)                           ![kB/baryon]
          ye(i) = state%bin%ye(iq)                                   ![]
          do it=1,nut
            ynu(it,i) = state%y(state%jy(iq,8+it))                   ![]
            tmp = (state%bin%n(iq)*state%y(state%jy(iq,10+it))**4)
     &        **(1./3.)
            if (state%y(state%jy(iq,10+it)).ge.0.) then
              znu(it,i) = tmp                              ![MeV/baryon]
            else
              znu(it,i) = -tmp                             ![MeV/baryon]
            endif
          enddo
        enddo

!.....isotropic diffusion source approximation..........................
        allocate(sout(n),edot(n),yout(n),ydot(n),pnu(n),vdot(n),
     &    ynudot(nut,n),znudot(nut,n))
        if (input_idsa%neutrino.eq.1) then
          iter = 20  !don't take more than 20 iterations
          iray = 1   !It is a 1-D call to nuprox
          call nuprox(n,dt,dV,d,s,ye,ynu,znu,iter,iray,
     &      sout,edot,yout,ydot,vdot,ynudot,znudot,dLmax)

!.....energy and pressure of trapped neutrinos..........................
!     note that contributions from the mu/tau neutrinos are neglected
          if (input_action%initial) call input_read_action
          pnu(1:n) = d*units%MeV/units%mb*(znu(len,:)+znu(lea,:))/3.
          tmp = 0.
          do i=1,n                                           ![g/cm/s^2]
            if (input_action%rel) then
              jrad(i) = 3.*pnu(i)/d(i)                       ![cm^2/s^2]
              do it=1,nut
                rate = src(:,it,i,iray)*d(i)/units%mb*dV(i)![particles/s]
                tmp = tmp + sum(rate*E)*units%MeV               ![erg/s]
              enddo
              hrad(i) = tmp/(4.*units%pi*ar(i)**2*d(i))      ![cm^3/s^3]
              krad(i) = pnu(i)                               ![g/cm/s^2]
            else
              jrad(i) = 0.
              hrad(i) = 0.
              krad(i) = 0.
            endif
          enddo
        else
          iter = 0
          edot = 0.
          ydot = 0.
          ynudot = 0.
          znudot = 0.
          dLmax = 0.
        endif

!.....substitute edot,ydot,vdot by parameterisation during collapse.....
        if (input_delept%parameterised.eq.1) then
          i = 1
          do while ((i.le.n).and.(d(i).gt.input_delept%d))
            if (s(i).gt.input_delept%s) then
              input_delept%parameterised = 0
              exit
            endif
            i = i+1
          enddo
        endif
        if (input_delept%parameterised.eq.1) then
          call delept_ysv(n,dt,dV,d,s,ye,edot,ydot,vdot)
        endif

!.....transfer lepton, energy, and momentum sources to hydrodynamics....
        vye = 0.
        venergy = 0.
        vstress = 0.
        do i=1,n-1
          iq = i+1
          vye(i) = ydot(i)*state%bda(iq)                          ![g/s]
          venergy(i) = edot(i)*state%bda(iq)                    ![erg/s]
          vstress(iq) = vdot(i)*state%ada(iq)                ![g*cm/s^2]
          do it=1,nut
            vynu(it,i) = ynudot(it,i)*state%bda(iq)               ![g/s]
            zold = state%y(state%jy(iq,10+it))                      ![?]
            tmp = znu(it,i) + znudot(it,i)*dt              ![MeV/baryon]
            if (tmp.ge.0.) then
              znew = (tmp**3/state%bin%n(iq))**0.25                 ![?]
            else
              znew = -((-tmp)**3/state%bin%n(iq))**0.25             ![?]
            endif
            vznu(it,i) = (znew-zold)/dt*state%bda(iq)           ![?*g/s]
          enddo
        enddo

!.....adjust time step to number of iterations..........................
        if (iter.gt.10) then
          dtnext = max(agile_red*dt,input_step%minimum)
        elseif (iter.le.5) then
          dtnext = min(agile_inc*dt,input_step%maximum)
        else
          dtnext = dt
        endif
        write(6,22) 'IDSA changes, iter: ',maxval(dLmax),iter
22      format(a24,2g12.4,i6)
        if (input_delept%parameterised.eq.1)
     &    write(6,*) 'Info: Parameterised deleptonisation'

        deallocate(ar,dV,d,s,ye,ynu,znu)
        deallocate(sout,edot,yout,ydot,pnu,vdot,ynudot,znudot)

        end subroutine agile_idsa

!=======================================================================

      end module agile_module

!***********************************************************************
