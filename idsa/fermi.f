!=======================================================================
!
!     IDSA: fermi, fermi distribution function
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      function fermi(arg)

      implicit none

      real, intent(in) :: arg
      real :: fermi

!-----------------------------------------------------------------------

      real :: tmp

      if (arg.gt.0.) then
        tmp = exp(-arg)
        fermi = tmp/(tmp+1.)
      else
        tmp = exp(arg)
        fermi = 1./(tmp+1.)
      endif

      end function fermi

!=======================================================================
