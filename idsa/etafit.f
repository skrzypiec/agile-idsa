!=======================================================================
!
!     IDSA: etafit, fitting Fermion temperature and degeneracy
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine etafit(nd,ed,beta,eta,status)

      use units_module

      implicit none

      integer, intent(out) :: status
      real, intent(in) :: nd,ed
      real, intent(out) :: beta,eta

!-----------------------------------------------------------------------
!
!     Input:
!     nd ... particle number density [1/cm^3]
!     ed ... particle energy density [MeV/cm^3]
!
!     Output:
!     beta ... inverse of kB*(particle temperature) [1/MeV]
!     eta ... particle degeneracy parameter []
!     status = 0 ... successful
!     status = 11 ... singular Jacobian
!     status = 12 ... negative eta in degenerate regime
!     status = 13 ... no convergence
!
!     This routine is based on the approximations of Fermi-Integrals
!     described in Epstein & Pethick, ApJ (1981).
!
!-----------------------------------------------------------------------

      logical, save :: initial=.true.
      integer, parameter :: itmax=20
      integer :: it
      real, parameter :: eps=1.d-06
      real, parameter :: fac=4.*units%pi/(units%h*units%c)**3
      real, parameter :: pi2=units%pi*units%pi, pi4=7.*pi2*pi2/120.
      real, parameter :: zeta=1.5*1.202	!3/2*zeta(3) Riemann zeta fct.
      real, parameter :: alpha=pi2/(3.*zeta)-1.
      real, parameter :: third=1./3.
      real, parameter :: etamin = -0.5*maxexponent(eps), etamax=1000.
      real :: nd3,ed4,eta2,zexp,F2,dF2,F3,dF3,tmp
      real, save :: F2max,F3max

!.....calculate Fermi integrals for etamax..............................
      if (initial) then
        eta2 = etamax*etamax
        zexp = zeta*exp(-alpha*etamax)
        F2max = (etamax*(eta2+pi2)/3. + zexp)**third
        zexp = pi4*exp(-etamax)
        F3max = (eta2*(eta2+2.*pi2)/4. + 2.*pi4 - zexp)**0.25
        initial = .false.
      endif
      it = 0	!initialize iteration counter

!.....check validity of number density..................................
      if (nd.le.0.) then	!negligible number density --> f=0
        status = 0
        goto 100
      endif

!.....check whether a fully degenerate solution is applicable...........
      nd3 = (nd/fac)**third
      tmp = nd3/F2max*F3max
      ed4 = (max(ed,0.)/fac)**0.25
      if (ed4.le.tmp) then      !energy density small --> degenerate limit
        eta = etamax
        beta = etamax/(3.**third*nd3)
        status = 0
        return
      endif

!.....try a non-degenerate solution assuming eta<0......................
      beta = pi4*nd/(zeta*ed)		!inverse temperature
      tmp = nd*beta**3/(fac*zeta)
      if (tmp.le.0.) then	!negligible number density --> f=0
        status = 0
        goto 100
      elseif (tmp.le.1.) then	!eta<0 consistent with assumption
        eta = log(tmp)
        status = 0
        return
      else		 	!take eta as guess for degenerate case
        eta = log(tmp)
      endif

!.....iterate to find degenerate solution assuming eta~>0...............
      do it=1,itmax
        eta2 = eta*eta
        zexp = zeta*exp(-alpha*eta)
        F2 = (eta*(eta2+pi2)/3. + zexp)**third
        dF2 = (eta2 + pi2/3. - alpha*zexp)/(3.*F2*F2)
        zexp = pi4*exp(-eta)
        F3 = (eta2*(eta2+2.*pi2)/4. + 2.*pi4 - zexp)**0.25
        dF3 = (eta*(eta2+pi2) + zexp)/(4.*F3*F3*F3)
        tmp = ed4*dF2 - nd3*dF3
        if (tmp.eq.0.) then	!singular Jacobian
          status = 11
          goto 100
        endif
        tmp = -(ed4*F2 - nd3*F3)/tmp
        eta = eta + tmp
!       write(6,*) 'eta,it: ',eta,it
        if (eta.lt.-0.5) then	!too far from assumption eta~>0
          status = 12
          goto 100
        elseif (abs(tmp).lt.eps*abs(eta)) then	!solution found
          eta2 = eta*eta
          zexp = zeta*exp(-alpha*eta)
          F2 = eta*(eta2+pi2)/3. + zexp
          beta = (fac*F2/nd)**third
          status = 0
          return
        endif
      enddo
      status = 13	!no convergence

!.....no solution found.................................................
100   continue
      if (status.ne.0) then
        write(6,*) 'Warning: no solution in etafit.f'
        write(6,*) 'nd,ed,status: ',nd,ed,status
      endif
      beta = 0.
      eta = etamin
      return

      end subroutine etafit

!=======================================================================
