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
      use ec_module, only: itmax,iymax,idmax,
     &                     ec_d,ec_t,ec_ye

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
      real :: Enu,tmp,etanp,etapn,egy,fe
      real :: t9
      real, dimension(ne) :: em_ec

!.TF.debug.eos-limits reached.
      if (t.le.0.5) then
        em(:,len)=0.
        ab(:,len)=0.
        em(:,lea)=0.
        ab(:,lea)=0.
        return
      endif
      if ((y(je).gt.0.5).or.(y(je).le.0.05)) then
        em(:,len)=0.
        ab(:,len)=0.
        em(:,lea)=0.
        ab(:,lea)=0.
        return
      endif

!.....nuclear ec-rates..................................................
      t9 = t/(units%kb/units%MeV)/1.e9
      if   (((t9.ge.ec_t(1)).and.(t9.le.ec_t(itmax))).and.
     &      ((y(je).ge.ec_ye(1)).and.(y(je).le.ec_ye(iymax))).and.
     &      ((d.ge.ec_d(1)).and.(d.le.ec_d(idmax))) )then
        call interpolate_ec(t9,y(je),d,em_ec(1:ne))
        em(1:ne,len)=em_ec(1:ne)*y(jh)*units%c
      else
        em(1:ne,len)=0.0
      endif

!.....nucleon degeneracy parameter......................................
      if (mu(jhat).gt.0.01) then
        etapn = (y(jn)-y(jp))/(exp(mu(jhat)/t)-1.0)*d/units%mb
        etanp = (y(jp)-y(jn))/(exp(-mu(jhat)/t)-1.0)*d/units%mb
      else
        etapn = y(jp)*d/units%mb
        etanp = y(jn)*d/units%mb
      endif

!.....energy dependence.................................................
      do ie=1,ne
        Enu = E(ie) + units%Q
        egy = Enu**2 * sqrt(1. - (units%me/Enu)**2)

!.....neutrino emissivity and absorptivity..............................
        fe=1.0/(1.0+exp((Enu-mu(je))/t))
        em(ie,len) = em(ie,len) + fac*etapn*egy*fe
        if (em(ie,len).le.0.0) then
          em(ie,len)=0.
          ab(ie,len)=0.
        else
          ab(ie,len)=em(ie,len)*(exp((Enu-mu(je)+mu(jhat))/t))
        endif

!.....threshold for positron capture....................................
        Enu = E(ie) - units%Q
        if (Enu-units%me.lt.0.) then
          em(ie,lea) = 0.
          ab(ie,lea) = 0.

!.....energy dependence.................................................
        else
          egy = Enu**2 * sqrt(1. - (units%me/Enu)**2)

!.....antineutrino emissivity and absorptivity..........................
          fe=1.0/(1.0+exp((Enu+mu(je))/t))
          em(ie,lea)=fac*etanp*egy*fe
          if (em(ie,lea).le.0.0) then
            em(ie,lea)=0.
            ab(ie,lea)=0.
          else
            ab(ie,lea)=em(ie,len)*(exp((Enu+mu(je)-mu(jhat))/t))
          endif
        endif
      enddo



      end subroutine nprates

!=======================================================================
