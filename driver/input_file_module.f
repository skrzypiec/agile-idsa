!***********************************************************************
!
!     DRIVER: input_file_module
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module input_file_module

      use agile_parameter_module

      implicit none

!=======================================================================
!
!     DRIVER: type definition for parameters in input file
!
!=======================================================================

      type input_action_type
        logical :: initial
        logical :: alloc      !flag for execution of allocation
        logical :: relax      !flag for execution of static relax
        integer :: step       !flag for time evolution and first step
        integer :: message    !flag that controls the amount of messages
        logical :: rel        !relativistic calculation
        logical :: restart    !restart flag (0...from progenitor, 1...from file no step)
      end type input_action_type

      type input_jacobi_type
        logical :: initial
        integer :: np        !number of coupled shells
        integer :: hoffset   !horizontal block offset
        integer :: voffset   !vertical row offset
        integer :: itmax     !maximal number of iterations
        real :: conv         !convergence
        real :: slowc        !damping for large errors
      end type input_jacobi_type

      type input_state_type
        logical :: initial
        integer :: np                     !number of progenitor zones
        integer :: nq                     !number of supernova zones
        integer :: ny                     !number of unknowns per zone
        character(24), dimension(nny) :: name  !name of variables
        integer, dimension(nny) :: scaling     !scaling mode
        real, dimension(nny) :: yscalmin       !minimum scale
        integer, dimension(nny) :: staggering  !staggering mode
      end type input_state_type
        
      type input_allocation_type
        logical :: initial
        logical :: grid                      !adaptive grid on/off
        real :: tau                          !time smoothing
        real :: alpha                        !space smoothing
        real :: visc			     !artificial viscosity
        real, dimension(nny) :: yloc         !influence on allocation
        integer, dimension(nny) :: scaling   !scaling mode
      end type input_allocation_type

      type input_step_type
        logical :: initial
        real :: eps			!precision time integration
        real :: start		        !start time
        real :: stop			!stop time
        real :: first			!size of first time step
        real :: minimum			!minimal acceptable time step
        real :: maximum			!step size limit
        real :: dt			!time interval for state_write
        integer :: dstep		!step interval for state_write
      end type input_step_type

      type input_idsa_type
        logical :: initial
        integer :: ne			!number of energy bins
	integer :: hydro		!do hydro steps
        integer :: neutrino		!do neutrino evolution
        integer :: mutau		!consider mu/tau neutrinos
      end type input_idsa_type          !during collapse phase

      type input_delept_type
        logical :: initial
        integer :: parameterised        !switch on or off
        real :: d                       !minimum density for entropy search
        real :: s                       !threshold entropy for off switch
      end type input_delept_type

      type input_path_type
        logical :: initial
        character(72) :: model          !path to initial profile and *.atb
        character(72) :: data           !path to input file and data output
        character(72) :: restart        !path to restart files
        character(72) :: eosfn          !file name for eos table
        character(72) :: progenitorfn	!file name for progenitor data
      end type input_path_type

!=======================================================================
!
!     DRIVER: parameters from the input file are stored here
!
!=======================================================================

      character(80) :: input_filename='filename'
      type(input_action_type) :: input_action
      type(input_jacobi_type) :: input_jacobi
      type(input_state_type) :: input_state
      type(input_allocation_type) :: input_allocation
      type(input_step_type) :: input_step
      type(input_idsa_type) :: input_idsa
      type(input_delept_type) :: input_delept
      type(input_path_type) :: input_path

      contains

!=======================================================================
!
!     DRIVER: input_initialisation
!
!=======================================================================

        subroutine input_initialisation

!.....set initial values here while waiting for f95.....................
        input_action%initial = .true.
        input_action%message = 1
        input_jacobi%initial = .true.
        input_state%initial = .true.
        input_allocation%initial = .true.
        input_step%initial = .true.
        input_step%start = 0.
        input_idsa%initial = .true.
        input_delept%initial = .true.
        input_path%initial = .true.

        end subroutine input_initialisation

!=======================================================================
!
!     DRIVER: input_read_action
!
!=======================================================================

        subroutine input_read_action

        integer :: ialloc,irelax,irel,irestart

        call input_search('ACTION    ')
        read(1,11) ialloc
        read(1,11) irelax
        read(1,11) input_action%step
        read(1,11) input_action%message
        read(1,11) irel
        read(1,11) irestart
11      format(t25,i12)
        close(1)
        input_action%alloc = (ialloc.eq.1)
        input_action%relax = (irelax.eq.1)
        input_action%rel = (irel.eq.1)
        input_action%restart = (irestart.eq.1)
        input_action%initial = .false.

        end subroutine input_read_action

!=======================================================================
!
!     DRIVER: input_read_jacobi
!
!=======================================================================

        subroutine input_read_jacobi

        call input_search('JACOBI      ')
        read(1,11) input_jacobi%np
        read(1,11) input_jacobi%hoffset
        read(1,11) input_jacobi%voffset
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,11) input_jacobi%itmax
        read(1,22) input_jacobi%conv
        read(1,22) input_jacobi%slowc
        close(1)
11      format(t25,i12)
22      format(t25,g12.4)
        input_jacobi%initial = .false.

        end subroutine input_read_jacobi

!=======================================================================
!
!     DRIVER: input_read_state
!
!=======================================================================
   
        subroutine input_read_state

        character(12) :: scaling,staggering
        integer :: iy

        call input_search('STATE     ')
        read(1,11) input_state%np
        read(1,11) input_state%nq
        read(1,11) input_state%ny
        read(1,*)
        read(1,*)
        read(1,*)
        do iy=1,input_state%ny
          read(1,22) input_state%name(iy),
     1      scaling,
     2      input_state%yscalmin(iy),
     3      staggering
          if (adjustl(scaling).eq.'logarithmic') then
            input_state%scaling(iy) = 2
          elseif (adjustl(scaling).eq.'sound_speed') then
            input_state%scaling(iy) = 3
          else
            input_state%scaling(iy) = 1
          endif
          if (adjustl(staggering).eq.'left') then
            input_state%staggering(iy) = -1
          elseif (adjustl(staggering).eq.'right') then
            input_state%staggering(iy) = 1
          else
            input_state%staggering(iy) = 0
          endif
        enddo
        close(1)
11      format(t25,i12)
22      format(a24,a12,g12.4,a12,g12.4)
        input_state%initial = .false.

        end subroutine input_read_state

!=======================================================================
!
!     DRIVER: input_read_allocation
!
!=======================================================================
   
        subroutine input_read_allocation

        character(12) :: grid,scaling
        integer :: iy

!.....first we need to know about the number of variables...............
        if (input_state%initial) call input_read_state

!.....then start reading the allocation data............................
        call input_search('ALLOCATION')
        read(1,11) grid
        input_allocation%grid = (adjustl(grid).eq.'on')
        if (input_allocation%grid) then
          write(6,*) 'Info: adaptive grid switched on'
        else
          write(6,*) 'Info: adaptive grid switched off'
        endif
        read(1,22) input_allocation%tau
        read(1,22) input_allocation%alpha
        read(1,22) input_allocation%visc
        read(1,*)
        read(1,*)
        read(1,*)
        do iy=1,input_state%ny
          read(1,22) input_allocation%yloc(iy),scaling
          if (adjustl(scaling).eq.'logarithmic') then
            input_allocation%scaling(iy) = 2
          else
            input_allocation%scaling(iy) = 1
          endif
        enddo
        close(1)
11      format(t25,a12)
22      format(t25,g12.4,a12,2g12.4)
        input_allocation%initial = .false.

        end subroutine input_read_allocation

!=======================================================================
!
!     DRIVER: input_read_step
!
!=======================================================================

        subroutine input_read_step

        call input_search('STEP      ')
        read(1,11) input_step%eps
        read(1,11) input_step%stop
        read(1,11) input_step%first
        read(1,11) input_step%minimum
        read(1,11) input_step%maximum
        read(1,11) input_step%dt
        read(1,22) input_step%dstep
        close(1)
11      format(t25,g12.4)
22      format(t25,i12)
        input_step%initial = .false.

        end subroutine input_read_step

!=======================================================================
!
!     IDSA: input_read_idsa
!
!=======================================================================

        subroutine input_read_idsa

        call input_search('IDSA      ')

        read(1,11) input_idsa%ne
        read(1,11) input_idsa%hydro
        read(1,11) input_idsa%neutrino
        read(1,11) input_idsa%mutau
        close(1)
11      format(t25,i12)
        input_idsa%initial = .false.
        
        end subroutine input_read_idsa

!=======================================================================
!
!     IDSA: input_read_delept
!
!=======================================================================

        subroutine input_read_delept

        call input_search('DELEPT    ')

        read(1,11) input_delept%parameterised
        read(1,22) input_delept%d
        read(1,22) input_delept%s
        close(1)
11      format(t25,i12)
22      format(t25,g12.4)
        input_delept%initial = .false.

        end subroutine input_read_delept

!=======================================================================
!
!     AGILE: input_read_path
!
!=======================================================================

        subroutine input_read_path

        call input_search('PATH      ')

        read(1,11) input_path%model
        read(1,11) input_path%data
        read(1,11) input_path%restart
        read(1,11) input_path%eosfn
        read(1,11) input_path%progenitorfn
11      format(t25,a72)
22      format(t25,i12)
        close(1)
        input_path%initial = .false.

        end subroutine input_read_path

!=======================================================================
!
!     DRIVER: input_search
!
!=======================================================================

        subroutine input_search(key)

        character(10), intent(in) :: key
        
!-----------------------------------------------------------------------
!
!     input:
!     key   search string
!
!     result:
!     the input file is opened and positioned at the first data entry
!     after the string key
!
!-----------------------------------------------------------------------

        character(10) :: word
        integer :: i

!.....open input-file...................................................
        open(1,file='model/'//trim(input_filename)//'.in',status='old')

!.....looking for keyword...............................................
        if (input_action%message.ge.2) write(6,*)
     1    'Info: reading input ',trim(key),' from ',
     2    trim(input_filename)//'.in'
        word = 'nokeyword '
        do while (key.ne.word)
          read(1,22,err=100) word
        enddo
22      format(a10)

!.....set positioner 4 rows further.....................................
        do i=1,4
          read(1,*)
        enddo

        return

!.....block not found...................................................
100     continue
        write(6,*) 'Error: block ',trim(key),' not found'
        write(6,*) 'in file ',trim(input_filename)//'.in!'
        stop

        end subroutine input_search

!=======================================================================

      end module input_file_module

!***********************************************************************
