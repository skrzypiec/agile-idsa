!=======================================================================
!
!     IDSA: rates between electron flavor neutrinos and free nucleons
!     Copyright: Andrey Yudin and Matthias Liebendoerfer,
!     30.09.2011, GNU public license.
!
!=======================================================================

      subroutine nprates(d,t,y,mu,em,ab)

      use units_module
      use species_module
      use egroup_module

      implicit none

      real, intent(in) :: d,t
      real, dimension(nuc+1), intent(in) :: y,mu
      real, dimension(ne,nut), intent(out) :: em,ab

!-----------------------------------------------------------------------
!
!     Andrey Yudin (2008) & Matthias Liebendoerfer (2004)
!
!     d ... density [g/cm^3]
!     t ... temperature [MeV]
!     y ... abundances [particles/baryon]
!     mu ... chemical potentials [MeV}
!     em ... neutrino emissivities [1/s]
!     ab ... neutrino absorptivities [1/s]
!
!-----------------------------------------------------------------------
 
!.....prefactors from units_module......................................
      real, parameter :: fac = (2.*units%pi/(units%h*units%c))**4
     &  * units%GF**2/units%pi * (units%gv**2+3.*units%ga**2) * units%c

      integer :: ie
      real :: Enu,tmp,etanp,etapn,egy,fexp

!.....special exponential factor (just to save computer time)...........
      fexp=exp(-mu(jhat)/t)

!.....nucleon degeneracy parameter......................................
      if (mu(jhat).gt.0.01) then
        etapn = (y(jn)-y(jp))*fexp/(1. - fexp)*d/units%mb
        etanp = (y(jp)-y(jn))/(fexp - 1.)*d/units%mb
      else
        etapn = y(jp)*d/units%mb
        etanp = y(jn)*d/units%mb
      endif

!.....energy dependence.................................................
      do ie=1,ne
        Enu = E(ie) + units%Q
        egy = Enu**2 * sqrt(1. - (units%me/Enu)**2)

!.....neutrino emissivity and absorptivity..............................
        tmp = (Enu-mu(je)+mu(jhat))/t	!units%Q is included in Enu
        tmp=exp(-tmp)
        ab(ie,len)=fac*etapn*egy/(tmp+fexp)
        em(ie,len)=ab(ie,len)*tmp

!.....threshold for positron capture....................................
        Enu = E(ie) - units%Q
        if (Enu-units%me.lt.0.) then
          em(ie,lea) = 0.
          ab(ie,lea) = 0.

!.....energy dependence.................................................
        else
          egy = Enu**2 * sqrt(1. - (units%me/Enu)**2)

!.....antineutrino emissivity and absorptivity..........................
          tmp = (Enu+mu(je)-mu(jhat))/t	!units%Q is included in Enu
          tmp=exp(-tmp)
          ab(ie,lea)=fac*etanp*egy*fexp/(fexp*tmp+1.)
          em(ie,lea)=ab(ie,lea)*tmp
        endif
      enddo

      end subroutine nprates

!=======================================================================
