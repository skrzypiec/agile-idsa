!***********************************************************************
!
!     IDSA: species module, particle properties
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module species_module

      implicit none

!.....particles in thermal equilibirum..................................
      integer, parameter :: nuc=4 !number of nuclear species
      integer, parameter :: jh=1  !heavy nucleus index
      integer, parameter :: ja=2  !alpha index
      integer, parameter :: jp=3  !proton index
      integer, parameter :: jhat=3 !if muhat is used instead of mup
      integer, parameter :: jn=4  !neutron index
      integer, parameter :: je=5  !electron-positron index

!.....size..............................................................
      real, dimension(nuc) :: A=(/0.,4.,1.,1./)
      real, dimension(nuc) :: Z=(/0.,2.,1.,0./)

!.....particles with energy spectrum....................................
      integer, parameter :: nut=2 !number of neutrinos
      integer, parameter :: len=1 !electron neutrino index
      integer, parameter :: lea=2 !electron antineutrino index

      end module species_module

!***********************************************************************
