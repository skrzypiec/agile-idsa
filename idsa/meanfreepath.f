!=======================================================================
!
!     IDSA: meanfreepath, calculates energy dependent mean free paths
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine meanfreepath(i,d,s,ye,ynu,znu,iray,t,f,feq,
     &  mfpe,mfpm,status)

      use egroup_module
      use eos_module
      use species_module
      use spectrum_module
      use units_module

      implicit none

      integer, intent(in) :: i,iray
      integer, intent(out) :: status
      real, intent(in) :: d,s,ye
      real, intent(out) :: t
      real, dimension(nut), intent(in) :: ynu,znu
      real, dimension(ne,nut), intent(out) :: f,feq,mfpe
      real, dimension(ne), intent(out) :: mfpm

!-----------------------------------------------------------------------
!
!     Input:
!     i ... zone identifier
!     d ... density [g/cm^3]
!     s ... specific entropy [kB/baryon]
!     ye ... electron fraction [particles/baryon]
!     ynu ... neutrino fractions [particles/baryon]
!     znu ... neutrino specific energy [MeV/baryon]
!     iray ... 1 for 1D call, 2 for 3D call
!
!     Output:
!     t                         ... matter temperature [MeV]
!     f(Energy,neutrino type)   ... neutrino distribution function
!                                   according to ynu and znu []
!     feq(Energy,neutrino type) ... equilibrium neutrino distribution
!                                   function []
!     mfpe(Energy,neutrino type) ... electron neutrino mean free path [cm]
!     mfpm(Energy)               ... mu/tau neutrino mean free path [cm]
!     status ... flag for error codes
!
!-----------------------------------------------------------------------

      integer :: ie,inu
      integer, dimension(nut) :: statnut
      real :: u
      real, dimension(nuc+1) :: y,mu,yscatt
      real, dimension(ne,nut) :: em,ab,chinp,chia

!.....construct energy groups...........................................
      if (egroup_empty) call egroup

!.....call equation of state............................................
      status = 0
      y(je) = ye
      call eosinterf(d,s,y,mu,t,u,status)
      if (status.ne.0) return

!.....determine neutrino distribution functions.........................
      call spectrum_add(i,d,ynu,znu,iray,f,statnut)
      do inu=1,nut
        if (statnut(inu).ne.0) then
          status = statnut(inu)
          return
        endif
      enddo

!.....calculate emission and absorption rates...........................
      call nprates(d,t,y,mu,em,ab)
      ab = em + ab	!stimulated absorption

!.....equilibrium distribution function.................................
      where (ab.le.0.)
        feq = 0.
      else where
        feq = em/ab
      end where

!.....calculate isoenergetic scattering kernels.........................
      call npscatt(d,t,y,chinp)
      call ionscatt(d,t,y,chia)

!.....calculate mean free path..........................................
      mfpe = units%c/(ab + chinp + chia)          !cm
      mfpm = units%c/(chinp(:,len) + chia(:,len)) !cm

      end subroutine meanfreepath

!=======================================================================
