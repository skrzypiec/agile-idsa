!=======================================================================
!
!     IDSA: uydot, calculates dfdt and sl
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine uydot(dt,d,t,y,mu,dr,ft,fs,etai,alpha,dfdt,sl)

      use egroup_module
      use species_module
      use spectrum_module
      use units_module

      implicit none

      real, intent(in) :: dt,d,t,dr
      real, dimension(nuc+1), intent(in) :: y,mu
      real, dimension(ne,nut), intent(in) :: ft,fs,etai,alpha
      real, dimension(ne,nut), intent(out) :: dfdt,sl

!-----------------------------------------------------------------------
!
!     dt ... time step [s]
!     d ... density [g/cm^3]
!     t ... temperature [MeV]
!     y ... abundances [particles/baryon]
!     mu ... chemical potentials [MeV]
!     dr ... propagation length [cm]
!     ft ... neutrino distribution functions of trapped neutrinos
!          ft = f(E) * 4.*pi/(h*c)**3*E**2*dE*mb/d  ![particles/baryon]
!     fs ... neutrino distribution functions of streaming neutrinos
!          fs = f(E) * 4.*pi/(h*c)**3*E**2*dE*mb/d  ![particles/baryon]
!     etai ... main factor in alpha for implicit finite differencing [1/s]
!     alpha ... diffusive neutrino source [particles/baryon/s]
!     dfdt ... rate of weak equilibrium adjustments [particles/baryon/s]
!     sl ... total approximate neutrino source [partcles/baryon/s]
!
!-----------------------------------------------------------------------
 
      real, parameter :: fac=4.*units%pi/(units%h*units%c)**3
      integer :: inu
      real, dimension(ne) :: ps
      real, dimension(ne,nut) :: em,ab,abr,sig

!.....calculate emission and absorption rates for transparent regime....
      call nprates(d,t,y,mu,em,ab)
      abr = ab  / (1. + ab*dt)
      abr = em + abr       !em*(1-f)-ab*f = em-(em+ab)*f            !1/s
      ab = em + ab

!.....phase space.......................................................
      ps = fac*E**2*dE*units%mb/d	![particles/baryon]
      do inu=1,nut
        em(:,inu) = ps*em(:,inu)		!particles/baryon/s
      enddo

!.....limit source function.............................................
      sig = ( (1.+ab*dt)*(alpha+abr*fs) + etai*dt*(em-ab*ft) )
     &  / (1. + (etai+ab)*dt)
      sig = min(max(sig,0.),em,ft*units%c/dr+abr*fs)
      dfdt = (em-ab*ft-sig)/(1. + ab*dt)
      sl = dfdt+sig-abr*fs
111   format(a12,20g12.4)

      end subroutine uydot

!=======================================================================
