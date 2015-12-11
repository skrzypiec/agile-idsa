!=======================================================================
!
!     IDSA: neutrino isoenergetic scattering rates
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine ionscatt(d,t,y,chi)

      use egroup_module
      use species_module
      use units_module

      implicit none

      real, intent(in) :: d,t
      real, dimension(nuc+1), intent(in) :: y
      real, dimension(ne,nut), intent(out) :: chi

!-----------------------------------------------------------------------
!
!     Input:
!     d ... density [g/cm^3]
!     t ... temperature [MeV]
!     y ... abundances [particles/baryon]
!
!     Output:
!     chi ... opacity [1/s]
!
!-----------------------------------------------------------------------

      real, parameter :: hvn=-0.5
      real, parameter :: hvp=0.5-2.0*units%sinsqthetw
      real, parameter :: c1=4.*units%pi**3*units%gf**2/units%h ![MeV*cm^6/s]
      real, parameter :: cv0=0.5*(hvp+hvn)
      real, parameter :: cv1=hvp-hvn
      real, parameter :: b=4.8d-06 ![1/MeV^2]
      real, parameter :: hc3 = (units%h*units%c)**3
      integer :: ie
      real :: na,nh,fa,fh,ya,yh,f0ya,f0yh,f1ya,f1yh,ris0,ris1

!.....ion densities.....................................................
      na = y(ja)*d/units%mb ![1/cm^3]
      nh = y(jh)*d/units%mb ![1/cm^3]

!.....coupling (note the minus sign)....................................
      fa = c1*na * (A(ja)*cv0 - (0.5*A(ja)-Z(ja))*cv1)**2 ![MeV*cm^3/s]
      fh = c1*nh * (A(jh)*cv0 - (0.5*A(jh)-Z(jh))*cv1)**2 ![MeV*cm^3/s]

!.....energy dependence.................................................
      do ie=1,ne

!.....alphas............................................................
        ya = 4.*b*A(ja)**(2./3.)*E(ie)**2	![]
        if (ya.gt.0.) then
          f0ya = (2.*ya - 1. + exp(-2.*ya))/ya**2
          f1ya = (2. - 3.*ya + 2.*ya**2 - (2.+ya)*exp(-2.*ya))/ya**3
        else
          f0ya = 0.
          f1ya = 0.
        endif

!.....nuclei............................................................
        yh = 4.*b*A(jh)**(2./3.)*E(ie)**2	![]
        if (yh.gt.0.) then
          f0yh = (2.*yh - 1. + exp(-2.*yh))/yh**2
          f1yh = (2. - 3.*yh + 2.*yh**2 - (2.+yh)*exp(-2.*yh))/yh**3
        else
          f0yh = 0.
          f1yh = 0.
        endif

!.....scattering kernel.................................................
        ris0 = fa*f0ya + fh*f0yh	![MeV*cm^3/s]
        ris1 = fa*f1ya + fh*f1yh	![MeV*cm^3/s]
        chi(ie,:) = 2.*(ris0 - ris1)/hc3*E(ie)**2 ![1/s]
      enddo

      end subroutine ionscatt

!=======================================================================
