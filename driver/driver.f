!=======================================================================
!
!     DRIVER AGILE-IDSA
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!   This program is free software: you can redistribute it and/or modify
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
!=======================================================================

      program driver

      use agile_module
      use eos_module, only: ls,eos_read
      use spectrum_module, only: spectrum_initialise
      use state_vector_module

      implicit none

!.....variables.........................................................
      character(len=72) agileinput
      character(len=4) :: fn
      logical :: writout
      integer :: i,count,status
      real :: dt,dthyd,dtnut,twrite
      type(state_vector) :: new

!.....initialise numerics...............................................
      call input_initialisation
      write(6,*) 'Name of input file? (without .in extension)'
      read(5,*) input_filename
      if (input_step%initial) call input_read_step
      if (input_action%initial) call input_read_action
      if (input_idsa%initial) call input_read_idsa
      if (input_delept%initial) call input_read_delept

!.....initialise equation of state and idsa.............................
      ls%initial = .true.
      if (input_path%initial) call input_read_path
      write(6,*) 'Info: reading '//trim(input_path%eosfn)
      call eos_read(input_path%eosfn)
      call spectrum_initialise

!-----build initial state-----------------------------------------------
      call agile_initialisation

!.....take time steps...................................................
      twrite = state%t          !time of last writeout
      istep = input_action%step !file number to start with
      dthyd = input_step%first  !first time step size
      dtnut = input_step%first  !first time step size
      count = 0                 !internal counter of time steps

!-----do time steps-----------------------------------------------------
      if (istep.gt.0) then
        do while (state%t.lt.input_step%stop)
          count = count + 1

!.....decide if data are written........................................
          if ((count.ge.input_step%dstep)
     &      .or.(state%t.ge.twrite+input_step%dt)) then
            writout = .true.
            count = 0
            twrite = state%t
          else
            writout = .false.
          endif

!.....do time step......................................................
          dt = min(dthyd,dtnut,input_step%stop-state%t)
          new = state
          status = 0
          call agile_step(dt,state,new,dthyd,equations_step)
          state = new
          call state_fill_zones(state,status)

!.....write hydro profile...............................................
          if (writout) then
            istep = istep+1
            call agile_write(istep)
          endif

!.....isotropic diffusion source approximation..........................
          call agile_idsa(dt,dtnut)

!.....prepare next time step............................................
          call state_scaling(state)
          if (input_action%message.ge.1) write(6,33) dt
33        format(1x,'Info: time step = ',g10.3)

!.....dump intermediate profiling information...........................
          if (input_action%message.ge.1) write(6,11)
11        format(72('-'))
        enddo

!.....write final result................................................
        istep = istep+1
        call agile_write(istep)
        write(6,22)
22      format(72('='))
      endif
      write(6,*) 'Info: done'

      end program driver

!=======================================================================
