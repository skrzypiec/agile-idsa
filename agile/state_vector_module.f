!***********************************************************************
!
!     AGILE: state_vector_module
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module state_vector_module

      use agile_parameter_module
      use input_file_module
      use thermo_module, only: thermo_input,thermo_output,thermo_eos
      use units_module

      implicit none

!=======================================================================
!
!     AGILE: type definition for state vector
!
!=======================================================================
!
!     The state vector contains the dimension nq of grid points and
!     the dimension ny of primitive variables. There is a big
!     1D array of length nq*ny that contains the data in a form
!     suitable for the Newton-Raphson scheme (which does not need
!     to distinguish the variables). The state vector also provides
!     separate arrays of length nq for the primitive and some 
!     derived variables. These are updated by the call to function
!     state_fill_zones and meant for easier access to the different
!     variables in the coding of the equations. The correspondence
!     between the two indexing convention is given by ic = jy(iq,iy),
!     where ic is the global index<=nq*ny and iq is the index of the
!     zone and iy the index of the variable.
!
!     Due to poor convergence when the zone widths are very small
!     with respect to the independent variable, the independent
!     variable is assigned to iy=1 and treated separately. A vector
!     dx(1) = x(1) and dx(iq) = x(iq)-x(iq-1) stores reference
!     value and increments of the independent variable. The value
!     of the independent variable is then given by
!     x(iq) + y(jy(iq,1)) as calculated by function state_sumdx.
!
!-----------------------------------------------------------------------

!.....state vector describing the whole state of the star...............
      type state_vector

!.....definition of unknowns............................................
        integer :: nq               	!number of grid points
        integer :: ny        		!# of unknowns per grid point
        integer, dimension(nnq,nny) :: jy 	!index transcription
        real, dimension(nnq) :: dx      !increments of independent var.
        real, dimension(nnc) :: y 	!unknowns
        real :: t			!time

!.....scaling...........................................................
        real, dimension(nnc) :: yscal   !individual scaling unknowns
        real, dimension(nny) :: yscl,ysclmin	!overall scales

!.....a copy of the state vector y......................................
        real, dimension(nnq) :: ax	!mass deviations [g]
        real, dimension(nnq) :: aa	!rest mass [g]
        real, dimension(nnq) :: ar 	!radius [cm]
        real, dimension(nnq) :: au 	!velocity [cm/s]
        real, dimension(nnq) :: am 	!gravitational mass [g]
        type(thermo_input) :: bin  	!input thermo quantities
        type(thermo_output) :: bout 	!output thermo quantities
        real, dimension(nnq) :: balpha	!lapse function []
        real, dimension(nnq,2) :: bynu  !neutrino fraction []
        real, dimension(nnq,2) :: bznu  !specific neutrino entropy [?]
					!(explained in header of state_nuput)

!.....further quantities................................................
        real, dimension(nnq) :: aV	!Volume [cm^3]
        real, dimension(nnq) :: aw	!4pi*r^2 * velocity/c [cm^2]
        real, dimension(nnq) :: amr	!m/r []
        real, dimension(nnq) :: aS	!"specific momentum" []
        real, dimension(nnq) :: ag	!gamma []
        real, dimension(nnq) :: bg	!gamma []
        real, dimension(nnq) :: aek	!"spec. kinetic energy"/c^2 []
        real, dimension(nnq) :: aep	!"spec. potential energy"/c^2 []
        real, dimension(nnq) :: bei	!"spec. internal energy"/c^2 []
        real, dimension(nnq) :: bp	!pressure [g/cm/s^2]
        real, dimension(nnq) :: bD	!"conserved density" [g/cm^3]
        real, dimension(nnq) :: ada	!mass differences on a [g]
        real, dimension(nnq) :: bda	!mass differences on b [g]
        real, dimension(nnq) :: brc     !radiation-matter coupling [1/s]
      end type state_vector

!.....overload assignments..............................................
      interface assignment (=)
        module procedure state_assignment
      end interface
        
      contains

!=======================================================================
!
!     AGILE: state_assignment
!
!=======================================================================

        subroutine state_assignment(new,old)

        type(state_vector), intent(in) :: old
        type(state_vector), intent(out) :: new

        integer :: n
        
        new%nq = old%nq
        new%ny = old%ny
        n = new%nq*new%ny
        new%jy(1:new%nq,1:new%ny) = old%jy(1:old%nq,1:old%ny)
        new%t = old%t
        new%dx(1:new%nq) = old%dx(1:old%nq)
        new%y(1:n) = old%y(1:n)

        end subroutine state_assignment

!=======================================================================
!
!     AGILE: state_scaling, updates the scales in the state vector
!
!=======================================================================

        subroutine state_scaling(state)

        type(state_vector) :: state

!-----------------------------------------------------------------------
!
!     input:
!     state   state vector
!
!     output:
!     state   state vector with updated scaling
!
!-----------------------------------------------------------------------

        integer :: n,ny,nq,iy
        integer, dimension(nnq) :: j
        real,save :: ln10
        real, dimension(nnq) :: x
        type(thermo_input) :: in
        type(thermo_output) :: out

!-----initialise--------------------------------------------------------

!.....read from input file..............................................
        if (input_state%initial) call input_read_state
        nq = state%nq
        ny = state%ny
        n = nq*ny
	ln10 = 1.d0/log(1.d1)

!-----relative scaling for odebs, newton and alloc----------------------

!.....independent variable..............................................
        iy = 1
        j(1:nq) = state%jy(1:nq,iy)

!.....default scaling...................................................
	if (input_state%scaling(iy).eq.1) then
          where(state%dx(1:nq).le.0.)
            state%yscal(j(1:nq)) = input_state%yscalmin(iy)
          else where
            state%yscal(j(1:nq)) = nq*state%dx(1:nq)
          end where

!.....logarithmic scaling...............................................
        elseif (input_state%scaling(iy).eq.2) then
	  state%yscal(j(1:nq)) = ln10
	else
          write(6,*) 'illegal scaling in state_scaling()!'
          stop
        endif

!.....loop over other variables.........................................
	do iy=2,ny
          j(1:nq) = state%jy(1:nq,iy)

!.....default scaling...................................................
	  if (input_state%scaling(iy).eq.1) then
	    state%yscal(j(1:nq)) = abs(state%y(j(1:nq))) 
     1        + input_state%yscalmin(iy)

!.....logarithmic scaling...............................................
          elseif (input_state%scaling(iy).eq.2) then
	    state%yscal(j(1:nq)) = ln10

!.....scaling of velocity by velocity of sound..........................
          elseif (input_state%scaling(iy).eq.3) then
            in%n(1:nq) = state%y(state%jy(1:nq,5))             ![g/cm^3]
            in%t(1:nq) = state%y(state%jy(1:nq,6))                ![MeV]
            in%ye(1:nq) = state%y(state%jy(1:nq,7))                  ![]
            call thermo_eos(in,out)
            where (out%status(1:nq).ne.0)
              out%cs(1:nq) = 0.
            end where
	    state%yscal(j(1:nq)) = abs(state%y(j(1:nq)))
     1        + input_state%yscalmin(iy) + out%cs(1:nq)          ![cm/s]
	  else
            write(6,*) 'illegal scaling in state_scaling()!'
            stop
          endif
        enddo

!-----scaling for deltaq------------------------------------------------

!.....set minimum scale.................................................
        state%ysclmin(1:ny) = input_state%yscalmin(1:ny)

!.....set range for variables...........................................
        call state_sumdx(state,x)
        state%yscl(1) = state%ysclmin(1) + sum(abs(x(2:nq)-x(1:nq-1)))
        do iy=2,ny
          j(1:nq) = state%jy(1:nq,iy)
          state%yscl(iy) = state%ysclmin(iy) + sum(abs(state%y(j(2:nq))
     1      - state%y(j(1:nq-1))))
        enddo
          
        end subroutine state_scaling

!=======================================================================
!
!     AGILE: state_fill_zones, updates physical quantities in zones
!
!=======================================================================

	subroutine state_fill_zones(s,status)

        integer, intent(out) :: status
        type(state_vector), intent(inout) :: s

!-----------------------------------------------------------------------
!
!     input:
!     s ... state vector
!
!     output:
!     status=0 ... zones updated successfully
!     s ... state vector with updated physical quantities
!
!-----------------------------------------------------------------------

        integer :: nq
        real, dimension(nnq) :: tmp

        status = 0
        nq = s%nq

!.....copy unknowns into physical variables.............................
        s%ax(1:nq) = s%y(s%jy(1:nq,1))       !mass deviations [g]
        call state_sumdx(s,s%aa)      !rest mass [g]
        s%ar(1:nq) = s%y(s%jy(1:nq,2))       !radius [cm]
        s%au(1:nq) = s%y(s%jy(1:nq,3))/units%c !velocity/c []
        s%am(1:nq) = s%y(s%jy(1:nq,4))       !grav mass [g]
        s%bin%n(1:nq) = s%y(s%jy(1:nq,5))    !rest mass density [g/cm^3]
        s%bin%t(1:nq) = s%y(s%jy(1:nq,6))    !temperature [MeV]
        s%bin%ye(1:nq) = s%y(s%jy(1:nq,7))   !electron fraction []
        s%balpha(1:nq) = s%y(s%jy(1:nq,8))   !lapse function []
        s%bynu(1:nq,1) = s%y(s%jy(1:nq,9))   !neutrino fraction []
        s%bynu(1:nq,2) = s%y(s%jy(1:nq,10))
        s%bznu(1:nq,1) = s%y(s%jy(1:nq,11))  !specific neutrino entropy [?]
        s%bznu(1:nq,2) = s%y(s%jy(1:nq,12))
        
!.....call equation of state............................................
        call thermo_eos(s%bin,s%bout)
        if (any(s%bout%status(1:nq).ne.0)) then
          status = 1
          return
        endif

!.....further quantities on the a-grid..................................
        tmp(1:nq) = 4.*units%pi*s%ar(1:nq)**2                    ![cm^2]
        s%aV(1:nq) = tmp(1:nq)*s%ar(1:nq)/3.                     ![cm^3]
        s%aw(1:nq) = tmp(1:nq)*s%au(1:nq)                        ![cm^2]
        s%amr(1) = 0.
        s%amr(2:nq) = units%G/units%c**2*s%am(2:nq)/s%ar(2:nq)       ![]
        if (input_action%rel) then
          tmp(1:nq) = 1. + s%au(1:nq)**2 - 2.*s%amr(1:nq)
          if (any(tmp(1:nq).le.0.)) then
            status = 1
            return
          endif
        else
          tmp(1:nq) = 1.
        endif
        s%ag(1:nq) = sqrt(tmp(1:nq))                                 ![]
        tmp(1:nq) = 2./(1.+s%ag(1:nq))
        s%aek(2:nq) = tmp(2:nq)*0.5*s%au(2:nq)**2                    ![]
        s%aek(1) = 0.1793*s%aek(2)
        s%aep(2:nq) = -tmp(2:nq)*s%amr(2:nq)                         ![]
        s%aep(1) = 0.1793*s%aep(2)
        
!.....further quantities on the b-grid..................................
        s%bp(1:nq) = s%bout%p(1:nq)                          ![g/cm/s^2]
        s%bda(2:nq) = s%dx(2:nq) + s%ax(2:nq) - s%ax(1:nq-1)        ![g]
        s%ada(1) = 0.5*s%bda(2)
        s%ada(2:nq-2) = 0.5*(s%bda(3:nq-1)+s%bda(2:nq-2))           ![g]
        s%ada(nq-1) = s%bda(nq) + 0.5*s%bda(nq-1)

!.....further quantities depending on both grids........................
        s%bg(1) = s%ag(1)
        s%bg(2:nq-1) = 0.5*(s%ag(2:nq-1) + s%ag(1:nq-2))             ![]
        s%bg(nq) = s%ag(nq)
        s%bD(1:nq) = s%bin%n(1:nq)/s%bg(1:nq)                  ![g/cm^3]
        s%bei(1:nq) = s%bout%e(1:nq)*s%bg(1:nq)                      ![]
        if (input_action%rel) then
          tmp(2:nq-1) = 0.5*(s%bout%e(3:nq)+s%bout%e(2:nq-1))
          tmp(nq) = s%bout%e(nq)
        else
          tmp(2:nq) = 0.
        endif
        s%aS(2:nq) = (1.+tmp(2:nq))*s%au(2:nq)                       ![]
        s%aS(1) = 0.2878*s%aS(2)

!.....radiation-matter coupling.........................................
        if (input_action%rel) then
          s%brc(1) = 0.
          s%brc(2:nq-1) = 4.*units%pi*units%G/units%c * s%ar(2:nq-1)
     1      * (s%bin%n(2:nq-1)*(1.+s%bout%e(2:nq-1))
     2      +  s%bout%p(3:nq)/units%c**2)                         ![1/s]
          s%brc(nq) = s%brc(nq-1)
        else
          s%brc = 0.
        endif

        end subroutine state_fill_zones

!=======================================================================
!
!       AGILE: update neutrino quantities
!
!=======================================================================

        subroutine state_nuput(ynu,znu,s)

        implicit none

        integer, parameter :: nut=2
        real, dimension(nut,nnq), intent(in) :: ynu,znu
        type(state_vector), intent(inout) :: s

!-----------------------------------------------------------------------
!
!       Input:
!       ynu(1) ... neutrino fraction [particle/baryon]
!       ynu(2) ... antineutrino fraction [particle/baryon]
!       znu(1) ... specific neutrino energy [MeV/baryon]
!       znu(2) ... specific antineutrino energy [MeV/baryon]
!
!       Output:
!       s ... statevector containing
!       y( 9)=bynu(1) ... neutrino fraction [particle/baryon]
!       y(10)=bynu(2) ... antineutrino fraction [particle/baryon]
!       y(11)=bznu(1) ... specific neutrino entropy [?]
!       y(12)=bznu(2) ... specific antineutrino entropy [?]
!
!       The units of the entropies are not relevant because in
!       equations_step only a conservation equation has to be solved.
!       In the Eulerian frame the conserved quantity is
!       (rho*y(11:12))**(3/4). However, for the Lagrangian scheme
!       in Agile we have to provide the specific conserved quantity, ie.
!       (y(11:12)**3/rho)**(1/4).
!
!-----------------------------------------------------------------------

        integer :: inu,iq
        real :: tmp

        do inu=1,nut
          do iq=1,s%nq
            tmp = ynu(inu,iq)
            s%y(s%jy(iq,8+inu)) = tmp
            s%bynu(iq,inu) = tmp
            if (znu(inu,iq).ge.0.) then
              tmp = (znu(inu,iq)**3/s%bin%n(iq))**0.25
            else
              tmp = -((-znu(inu,iq))**3/s%bin%n(iq))**0.25
            endif
            s%y(s%jy(iq,10+inu)) = tmp
            s%bznu(iq,inu) = tmp
          enddo
        enddo

        end subroutine state_nuput

!=======================================================================
!
!       AGILE: update neutrino quantities
!
!=======================================================================

        subroutine state_nuget(s,ynu,znu)

        implicit none

        integer, parameter :: nut=2
        type(state_vector), intent(in) :: s
        real, dimension(nut,nnq), intent(out) :: ynu,znu

!-----------------------------------------------------------------------
!
!       Input:
!       s ... statevector containing
!       y( 9)=bynu(1) ... neutrino fraction [particle/baryon]
!       y(10)=bynu(2) ... antineutrino fraction [particle/baryon]
!       y(11)=bznu(1) ... specific neutrino entropy [?]
!       y(12)=bznu(2) ... specific antineutrino entropy [?]
!
!       Output:
!       ynu(1) ... neutrino fraction [particle/baryon]
!       ynu(2) ... antineutrino fraction [particle/baryon]
!       znu(1) ... specific neutrino energy [MeV/baryon]
!       znu(2) ... specific antineutrino energy [MeV/baryon]
!
!       The units of the entropies are not relevant because in
!       equations_step only a conservation equation has to be solved.
!       In the Eulerian frame the conserved quantity is
!       (rho*y(11:12))**(3/4). However, for the Lagrangian scheme
!       in Agile we have to provide the specific conserved quantity, ie.
!       (y(11:12)**3/rho)**(1/4).
!
!-----------------------------------------------------------------------

        real, parameter :: third=1./3.
        integer :: inu,iq
        real :: tmp

        do inu=1,nut
          do iq=1,s%nq
            ynu(inu,iq) = s%y(s%jy(iq,8+inu))
            tmp = (s%bin%n(iq)*s%y(s%jy(iq,10+inu))**4)**third
            if (s%y(s%jy(iq,10+inu)).ge.0.) then
              znu(inu,iq) = tmp
            else
              znu(inu,iq) = -tmp
            endif
          enddo
        enddo

        end subroutine state_nuget

!=======================================================================
!
!     AGILE: read state
!
!=======================================================================

        subroutine state_read(kount,new)
        
        integer, intent(in) :: kount
        type(state_vector), intent(out) :: new
        
!-----------------------------------------------------------------------
!
!     input:
!     kount   number in filename
!
!     output:
!     new   state read from file
!
!-----------------------------------------------------------------------

        character(58) :: name
        character(24), dimension(nny) :: comment
        character(4) :: ext,num
        integer :: iy,iq,i,iy0,iy1,nq
        real :: dummy
        real, dimension(nnq) :: x

!.....initialise........................................................
        if (input_action%initial) call input_read_action
      
!.....generate filename.................................................
        if (kount.lt.-2) then
          stop 'unallowed kount in getstep.f!'
        elseif (kount.eq.-2) then
          ext = '.prf'
          num = '   '
        elseif (kount.eq.-1) then
          ext = '.dmp'
          num = '   '
        else
          ext = '.stp'
          num = char(48+mod(kount/1000,10))
     1      //char(48+mod(kount/100,10))
     2      //char(48+mod(kount/10,10))
     3      //char(48+mod(kount,10))
        endif
        if (input_path%initial) call input_read_path
        name = trim(input_path%data)
     1       //trim(input_filename)//trim(num)//ext
        if (input_action%message.ge.1) write(6,11) trim(name)
11      format(1x,'Info: reading ',a)

!.....read header.......................................................
        open(1,file=name,status='old')
        do i=1,2
          read(1,*)
        enddo
        read(1,22) new%nq,new%ny
22      format(t16,2i15)
        new%ny = new%ny-1
        read(1,33) new%t
33      format(t16,g24.16)
        do i=1,3
          read(1,*)
        enddo
        read(1,44) (comment(iy)(1:12), iy=1,new%ny)
        read(1,44) (comment(iy)(13:24), iy=1,new%ny)
44      format(t17,1000(a12,3x))
        do i=1,2
          read(1,*)
        enddo

!.....read data.........................................................
        nq = new%nq
        iy0=1
        iy1 = new%ny
        do iq=1,nq
          new%jy(iq,1:new%ny) = (/(i,i=iy0,iy1)/)
          read(1,55) new%y(iy0:iy1)
          iy0 = iy1+1
          iy1 = iy1 + new%ny
        enddo
55      format(t16,1000g15.7)

!.....read dx..........................................................
        do i=1,3
          read(1,*)
        enddo
        do iq=1,nq
          read(1,99) dummy,new%dx(iq),new%y(new%jy(iq,1))
        enddo
99      format(g15.7,2g30.22)
        close(1)
      
        end subroutine state_read

!=======================================================================
!
!     AGILE: write state
!
!=======================================================================

        subroutine state_write(kount,new)

        integer, intent(in) :: kount
        type(state_vector), intent(inout) :: new
        
!-----------------------------------------------------------------------
!
!     input:
!     kount   number in filename
!     new     state vector
!
!-----------------------------------------------------------------------

        character(58) :: name
        character(24) :: comment
        character(4) :: ext,num
        integer :: iy,iq,n
        real, dimension(nnq) :: x

!.....initialise........................................................
        if (input_action%initial) call input_read_action
        if (input_state%initial) call input_read_state      

!.....prepare data......................................................
        n = new%nq*new%ny
        where (abs(new%y(1:n)).lt.1.e-99)
          new%y(1:n) = 0.
        end where

!.....generate filename.................................................
        if (kount.lt.-2) then
          stop 'unallowed kount in writstep.f!'
        elseif (kount.eq.-2) then
          ext = '.prf'
          num = '   '
        elseif (kount.eq.-1) then
          ext = '.dmp'
          num = '   '
        elseif (kount.lt.10000) then
          ext = '.stp'
          num = char(48+mod(kount/1000,10))
     1      //char(48+mod(kount/100,10))
     2      //char(48+mod(kount/10,10))
     3      //char(48+mod(kount,10))
        else
          write(6,*) 'Error: 10000 dumps reached in writstep.f!'
          stop
        endif
        if (input_path%initial) call input_read_path
        name = trim(input_path%data)
     1       //trim(input_filename)//trim(num)//ext
        if (input_action%message.ge.1) then
          write(6,99) new%t
99        format(1x,'Info: time = ',g24.16)
          write(6,11) name
11        format(1x,'Info: writing ',a58)
        endif

!.....write header......................................................
        open(1,file=name,status='unknown')
        do iy=1,new%ny
          write(1,22,advance='no')
        enddo
        write(1,22)
22      format('===============')
        write(1,33) adjustl(name)
33      format(1x,'file:',t16,a58)
        write(1,44) new%nq,new%ny+1
44      format(1x,'size:',t16,2i15)
        write(1,55) new%t
55      format(1x,'time:',t16,g24.16)
        do iy=1,new%ny
          write(1,22,advance='no')
        enddo
        write(1,22)
        write(1,*)
        do iy=1,new%ny
          write(1,66,advance='no')
        enddo
        write(1,66)
66      format('|--------------')
        write(1,77) 'shell       ',
     1    (input_state%name(iy)(1:12),iy=1,new%ny)
77      format(1x,1000(a12,3x))
        write(1,77) 'number      ',
     1    (input_state%name(iy)(13:24),iy=1,new%ny)
        do iy=1,new%ny
          write(1,66,advance='no')
        enddo
        write(1,66)
        write(1,*)

!.....write data........................................................
        call state_sumdx(new,x)
        do iq=1,new%nq
          write(1,88) float(iq),x(iq),new%y(new%jy(iq,2:new%ny))
        enddo
88      format(1000e15.7)
        write(1,*)
        do iy=1,new%ny
          write(1,22,advance='no')
        enddo
        write(1,22)

!.....write dx.........................................................
        write(1,*)
        do iq=1,new%nq
          write(1,111) float(iq),new%dx(iq),new%y(new%jy(iq,1)),
     &      new%bout%s(iq),new%bout%cs(iq),new%bout%status(iq)
        enddo
111     format(e15.7,2e30.22,2e15.7,i5)
        close(1)

        end subroutine state_write

!=======================================================================
!
!     state_sumdx
!
!=======================================================================

        subroutine state_sumdx(state,x)

        implicit none

        real, dimension(nnq), intent(out) :: x
        type(state_vector), intent(in) :: state

!-----------------------------------------------------------------------
!
!     Calculates the independent variable from the increments stored
!     in state%dx. A flag chooses whether to include y(jy(:,1)).
!
!     Input:
!     state ... state vector
!
!     Output:
!     x ... sum(state%dx)+y(jy(:,1))
!
!-----------------------------------------------------------------------

        integer :: iq,nq

        nq = state%nq
        x(1) = state%dx(1)
        do iq=2,nq
          x(iq) = x(iq-1) + state%dx(iq)
        enddo
        x(1:nq) = x(1:nq) + state%y(state%jy(1:nq,1))

        end subroutine state_sumdx

!=======================================================================
!
!     state_newdx
!
!=======================================================================

        subroutine state_newdx(state)

        implicit none

        type(state_vector), intent(inout) :: state

!-----------------------------------------------------------------------

        integer :: iq,nq
        real, dimension(nnq) :: y

        nq = state%nq
        y(1:nq) = state%y(state%jy(1:nq,1))
        state%dx(1) = state%dx(1) + y(1)
        do iq=2,nq
          state%dx(iq) = state%dx(iq) + (y(iq)-y(iq-1))
        enddo
        state%y(state%jy(1:nq,1)) = 0.
        do iq=2,nq
          if (state%dx(iq).lt.0.) then
            write(6,*) 'Error: dx<0 in zone ',iq
            stop
          endif
        enddo

        end subroutine state_newdx

!=======================================================================

      end module state_vector_module

!***********************************************************************
