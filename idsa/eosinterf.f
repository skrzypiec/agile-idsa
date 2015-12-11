!=======================================================================
!
!     IDSA-EOS interface
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine eosinterf(d,s,y,mu,t,u,status)

      use eos_module
      use species_module
      use units_module

      implicit none

      integer :: status
      real, intent(in) :: d
      real, intent(in) :: s
      real, dimension(nuc+1) :: y,mu
      real, intent(out) :: t,u

!-----------------------------------------------------------------------
!
!     ML, CITA (2004)
!     d ... density [g/cm^3]
!     s ... entropy [kB/baryon]
!     y ... abundances [particles/baryon]
!     mu ... chemical potentials [MeV]
!     t ... temperature [MeV]
!     u ... specific interal energy [MeV/baryon]
!
!-----------------------------------------------------------------------

      real :: cs
      type(eos_state) :: state

!.....call equation of state............................................
      call eos_sinterp(d,y(je),s,cs,state,status)
      if (status.ne.0) return

!.....eliminate negative abundances.....................................
      state%xa = max(0.,state%xa)
      state%xp = max(0.,state%xp)
      state%xn = max(0.,state%xn)

!.....set output variables..............................................
      A(jh) = state%a
      Z(jh) = state%z
      if (A(jh).gt.0.) then
        y(jh) = max(1.-state%xa-state%xp-state%xn,0.)/A(jh) !heavy nuclei
      else
        y(jh) = 0.
      endif
      y(ja) = state%xa/A(ja) !alpha abundance
      y(jp) = state%xp/A(jp) !proton abundance
      y(jn) = state%xn/A(jn) !neutron abundance
      mu(jhat) = state%muhat !muhat without mass difference
      mu(jn) = state%mun !not used
      mu(je) = state%mue !electron chemical potential
      t = state%t !temperature
      u = state%e*units%mb/units%MeV !specific internal energy

!.....check whether all abundances are positive.........................
      if (any(y.lt.0.)) then
        write(6,*) 'Warning: negative abundances in eosinterf.f'
        write(6,11) 'input: ',d,y(je),s
        write(6,11) 'output: ',state
11      format(a8,12g12.4)
        status = 5
      endif
      
      end subroutine eosinterf

!=======================================================================
