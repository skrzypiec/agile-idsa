!=======================================================================
!
!     IDSA: neutrino isoenergetic scattering rates
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine npscatt(d,t,y,chi)

      use units_module
      use species_module
      use egroup_module

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

      real, parameter :: hbar = 0.5*units%MeV*units%h/units%pi ![erg*s]
      real, parameter :: c1 = 0.5*hbar*hbar/units%mb ![erg*cm^2]
      real, parameter :: hvn=-0.5
      real, parameter :: han=-0.5*units%ga
      real, parameter :: hvp=0.5-2.0*units%sinsqthetw
      real, parameter :: hap=0.5*units%ga
      real, parameter :: c2=8.*units%pi**3*units%gf**2/units%h ![MeV*cm^6/s]
      real, parameter :: cn0=hvn*hvn + 3.*han*han !Legendre coeff 0
      real, parameter :: cp0=hvp*hvp + 3.*hap*hap
      real, parameter :: cn1=(hvn*hvn - han*han)/3. !Legendre coeff 1
      real, parameter :: cp1=(hvp*hvp - hap*hap)/3.
      real, parameter :: hc3 = (units%h*units%c)**3
      integer :: ie
      real :: nn,np,efern,eferp,tmp,etann,etapp,ris0,ris1

!.....nucleon particle densities........................................
      nn = y(jn)*d/units%mb ![1/cm^3]
      np = y(jp)*d/units%mb ![1/cm^3]

!.....degeneracy parameter eta..........................................
      if (nn.le.0.) then
        etann = 0.
      else
        efern = c1*(3.*units%pi**2*nn)**(2./3.)/units%MeV ![MeV]
        tmp = 1.5*t/efern ![]
        etann = nn*tmp/sqrt(1.+tmp**2) ![1/cm^3]
      endif
      if (np.le.0.) then
        etapp = 0.
      else
        eferp = c1*(3.*units%pi**2*np)**(2./3.)/units%MeV ![MeV]
        tmp = 1.5*t/eferp ![]
        etapp = np*tmp/sqrt(1.+tmp**2) ![1/cm^3]
      endif

!.....energy dependence.................................................
      do ie=1,ne
        ris0 = c2*(etann*cn0 + etapp*cp0) ![MeV*cm^3/s]
        ris1 = c2*(etann*cn1 + etapp*cp1) ![MeV*cm^3/s]
        chi(ie,:) = 2.*(ris0 - ris1)/hc3*E(ie)**2 ![1/s]
      enddo

      end subroutine npscatt

!=======================================================================
