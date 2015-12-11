!***********************************************************************
!
!     IDSA: egroup_module, define neutrino energy groups
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module egroup_module

      implicit none

      logical :: egroup_empty=.true.
      integer, parameter :: ne=20
      real, dimension(ne) :: e,de

      contains

!=======================================================================
!
!     IDSA: egroup
!
!=======================================================================

      subroutine egroup

!-----------------------------------------------------------------------
!
!     Geometric series with E^2dE = d(E^3) as suggested by S. W. Bruenn
!
!-----------------------------------------------------------------------

      real, parameter :: emin=3.,emax=300.

      integer :: ie
      real :: f1,f2

      f1 = (emax/emin)**(1./float(ne-1))
      f2 = (f1-1.)/sqrt( (1.+f1*(1.+f1))/3. )
      e(1) = emin
      de(1) = f2*e(1)
      do ie=2,ne
        e(ie) = f1*e(ie-1)
        de(ie) = f2*e(ie)
      enddo
      egroup_empty = .false.

      end subroutine egroup

!=======================================================================

      end module egroup_module

!***********************************************************************
