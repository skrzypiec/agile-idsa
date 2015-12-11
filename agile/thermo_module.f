!***********************************************************************
!
!     AGILE-EOS interface
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module thermo_module

      use agile_parameter_module
      use eos_module, only: eos_state, eos_tinterp
      use input_file_module
      use units_module

      implicit none

!=======================================================================
!
!     Type definition for AGILE thermo quantities
!
!=======================================================================

      type thermo_input
        real, dimension(nnq) :: n   !rest mass density [g/cm^3]
        real, dimension(nnq) :: t   !temperature [MeV]
        real, dimension(nnq) :: ye  !electron fraction []
      end type thermo_input

      type thermo_output
        integer, dimension(nnq) :: status !status flag
        real, dimension(nnq) :: p   !pressure [g/cm/s^2]
        real, dimension(nnq) :: e   !(specific energy)/c^2 []
        real, dimension(nnq) :: s   !entropy per baryon [kB/baryon]
        real, dimension(nnq) :: cs  !sound speed [cm/s]
      end type thermo_output

      contains

!=======================================================================
!
!     AGILE-EOS interface
!
!=======================================================================

        subroutine thermo_eos(in,out)

        type(thermo_input), intent(in) :: in
        type(thermo_output), intent(out) :: out
        
!-----------------------------------------------------------------------
!
!     Input:
!     n ... rest mass density [g/cm^3]
!     temp ... temperature [MeV]
!     ye ... electron fraction []
!
!     Output:
!     status ... status flag reporting on success of eos call
!     p ... pressure [g/cm/s^2]
!     e ... (specific energy)/c^2 []
!     s ... entropy [kB/baryon]
!     cs ... velocity of sound [cm/s]
!
!-----------------------------------------------------------------------

        integer :: iq,status
        real :: d,t,ye,cs
        type(eos_state) :: state

!.....loop over zones...................................................
        if (input_action%initial) call input_read_action
        if (input_state%initial) call input_read_state
        do iq=1,input_state%nq

!.....convert AGILE input quantities to cgs.............................
          d = in%n(iq)                                      ![g/cm^3]
          t = in%t(iq)                                      ![MeV]
          ye = in%ye(iq)                                    ![]

!.....call equation of state............................................
          call eos_tinterp(d,ye,t,cs,state,status)
          if (status.ne.0) then
            write(6,*) 'Warning: eos status = ',status,' in zone ',iq
            write(6,*) 'd/t/ye: ',d,t,ye
          endif

!.....convert output quantities to AGILE units..........................
          out%p(iq) = state%p                               ![g/cm/s^2]
          out%e(iq) = state%e/units%c**2                    ![]
          out%s(iq) = state%s                               ![kB/baryon]
          out%cs(iq) = cs                                   ![cm/s]
        enddo

        end subroutine thermo_eos

!=======================================================================

      end module thermo_module

!***********************************************************************
