!=======================================================================
!
!     AGILE: initialisation, i.e. grid point allocation and relax
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine agile_start(new,equations_relax)

      use adaptive_grid_module, only: adaptive_linint_template
      use input_file_module
      use newton_module
      use state_vector_module
      
      implicit none

      external equations_relax

      type(state_vector), intent(out) :: new

!-----------------------------------------------------------------------
!
!     input:
!     equations_relax   routine with equations for static relax
!
!     output: 
!     new
!
!     if (alloc): make a grid allocation guess --> write *0000.stp
!     if (relax): read *0000.stp --> relax grid --> write *0001.stp
!     if (step): prepare state vector for time steps
!
!-----------------------------------------------------------------------

      integer :: status=0
      integer :: itmax
      real :: slowc,conv
      type(state_vector) :: old
      
!.....initialise........................................................
      if (input_action%initial) call input_read_action

!.....restart from step.................................................
      if (input_action%restart) then
        write(6,*) input_action%step
        call state_read(input_action%step,old)
        new = old

      else

!.....read and rewrite stellar profile..................................
      call progenitor(old)
      call state_write(-2,old)
      adaptive_linint_template = old
      new = old

!-----grid point allocation---------------------------------------------
      if (input_action%alloc) then

!.....allocate grid points (if adaptive grid is switched on)............
        if (input_allocation%initial) call input_read_allocation
        if (input_allocation%grid) then
          call grid_initialisation(old,new,status)
          if (status.ne.0) then
            write(6,*) 'Error: grid point allocation failed'
            call state_write(-1,old)
            stop
          endif
        endif
        
!.....write state with allocated grid points............................
        if (input_step%initial) call input_read_step
        new%t = input_step%start
        call state_write(0,new)
        write(6,11)
11      format(72('='))
      endif
      
!-----relax-------------------------------------------------------------
      if (input_action%relax) then

!.....read in allocation guess..........................................
        call state_read(0,old)

!.....complete old state vector.........................................
        call state_fill_zones(old,status)
        if (status.ne.0) return
        call state_scaling(old)

!.....relaxation........................................................
        itmax = 20
        slowc = 0.1d+0
        conv = 0.1d-8
        newton_flag%block = .true.
        newton_flag%analyze = .true.
        newton_flag%banddiag = .true.
        newton_flag%bigrad = .false.
        new = old
        call newton_raphson(itmax,slowc,conv,old,
     1    new,status,equations_relax)
        if (status.ne.0) then
          write(6,*) 'Error: relax failed'
          call state_write(-1,old)
          stop
        else
          if (input_action%message.ge.1)
     1      write(6,*) 'Info: successful relax'
        endif
      
!.....write relaxed state...............................................
        call state_write(1,new)
        write(6,11)
      endif

      endif

!-----fill state vector-------------------------------------------------
      if (input_action%step.gt.0) then
        call state_read(input_action%step,new)
        status = 0
        call state_fill_zones(new,status)
        if (status.ne.0) then
          write(6,*) 'Error: invalid zone values in initial state!'
          call state_write(-1,new)
          stop
        endif
        call state_scaling(new)
      endif

      
      end subroutine agile_start

!=======================================================================
