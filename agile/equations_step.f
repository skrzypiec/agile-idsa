!=======================================================================
!
!     AGILE: Equations for hydro time step
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine equations_step(jacobi,o,n,f,status)

      use agile_module, only: vstress,venergy,vye,vynu,vznu,
     &  jrad,hrad,krad
      use agile_parameter_module
      use state_vector_module

      implicit none

      logical, intent(in) :: jacobi
      integer, intent(out) :: status
      real, dimension(nnc), intent(out) :: f
      type(state_vector), intent(in) :: o,n

!-----------------------------------------------------------------------
!
!     input:
!     jacobi   is true if the routine is called for derivatives
!     o   old state vector
!     n   guess state vector
!
!     output:
!     f   residua
!     status is set to 1 in case of a problem (not used)
!
!     This subroutine implements the hydrodynamics equations according
!     to Eqs. (88)-(95) in Liebendoerfer et al., ApJS 150 (2004).
!
!     Notes:
!     a) The third line of Eq. (90) in the paper contains a typo:
!     wrong: (2u^2-m/r)*2p/rho, correct: (2u^2-m/r)*p/rho.
!     b) venergy is defined here as MINUS the integral shown
!     in the paper, i.e. venergy=-e^{ext} in Eq. (96).
!     c) The terms handling the viscosity Q are not included in the
!     equations of the 2004 paper. They are described in Liebendoerfer,
!     Rosswog & Thielemann, ApJS 141 (2002).
!
!-----------------------------------------------------------------------

      logical, save :: initial=.true.
      logical, dimension(nnq) :: interior
      integer, save :: nq,nq1,nz,nz1
      integer :: ie,i,inu
      integer, dimension(nnq+1,nny), save :: jf
      real :: dt,opsurf
      real, save :: psurf
      real, dimension(nnq) :: aurel,burel,divu,comp,bQ,baVQ,bdr
      real, dimension(nnq) :: atmp,btmp,grad,source
      real, dimension(nnq) :: dSdt,dDdt,deidt,dEdt,dyedt,advect,dnudt
      real, dimension(nnq) :: Sflux,Dflux,aEflux,bEflux,yeflux
      real, dimension(nnq,4) :: ynuflux,znuflux
      real, dimension(nnq) :: bdV,bdw,bp,aalpha,aQ,bdu,bV,apq,an,ae,eta
      real, dimension(nnq) :: avstress,bvenergy,bvye,obdV
      real, dimension(nnq) :: bjrad,bhrad,bkrad,bradcop,arad

!-----initialise--------------------------------------------------------

!.....abbreviations for array-lengths...................................
      if (initial) then
        nq = n%nq
        nq1 = nq-1
        nz = nq+1
        nz1 = nz-1

!.....assign zones to residuum vector...................................
        if (input_jacobi%initial) call input_read_jacobi
        do ie=input_jacobi%voffset+1,n%ny
          jf(1:nz1,ie) = (/(i,i=0,nz-2)/)*n%ny+ie-input_jacobi%voffset
        enddo
        do ie=1,input_jacobi%voffset
          jf(2:nz,ie) = (/(i,i=1,nz1)/)*n%ny+ie-input_jacobi%voffset
        enddo
      
!.....fixed surface pressure............................................
        psurf = o%bout%p(nq)                                 ![g/cm/s^2]
        if (input_idsa%initial) call input_read_idsa
        if (input_allocation%initial) call input_read_allocation
        initial = .false.
      endif

!-----precalculate some terms for later use-----------------------------
      status = 0
      dt = n%t - o%t                                                ![s]
      bV(1) = n%aV(1)
      bV(2:nq1) = 0.5*(n%aV(2:nq1) + n%aV(1:nq-2))               ![cm^3]
      bV(nq) = n%aV(nq)
      bdr(2:nq) = n%ar(2:nq) - n%ar(1:nq1)                         ![cm]
      bdV(2:nq) = 4./3.*units%pi*bdr(2:nq)                       ![cm^3]
     1  * (n%ar(2:nq)**2 + n%ar(1:nq1)**2 + n%ar(2:nq)*n%ar(1:nq1))
      obdV(2:nq) = 4./3.*units%pi*(o%ar(2:nq) - o%ar(1:nq1))     ![cm^3]
     1  * (o%ar(2:nq)**2 + o%ar(1:nq1)**2 + o%ar(2:nq)*o%ar(1:nq1))
      bdu(2:nq) = n%au(2:nq) - n%au(1:nq1)                           ![]
      bdw(2:nq) = n%aw(2:nq) - n%aw(1:nq1)                       ![cm^2]
      bp(1:nq1) = n%bp(1:nq1)                                ![g/cm/s^2]
      aalpha(1) = n%balpha(1)
      aalpha(2:nq1) = 0.5*(n%balpha(3:nq) + n%balpha(2:nq1))         ![]
      aalpha(nq) = n%balpha(nq)
    
!.....set surface pressure..............................................
      bp(nq) = psurf                                         ![g/cm/s^2]
      opsurf = psurf                                         ![g/cm/s^2]

!.....convert neutrino influence........................................
      bjrad(2:nq1) = jrad(1:nq-2)/units%c**2                         ![]
      bhrad(2:nq1) = hrad(1:nq-2)/units%c**3                         ![]
      bkrad(2:nq1) = krad(1:nq-2)                            ![g/cm/s^2]
      bjrad(nq) = 0.
      bhrad(nq) = 0.
      bkrad(nq) = 0.
      avstress(2:nq1) = aalpha(2:nq1)*vstress(2:nq1)/units%c      ![g/s]
      bvenergy(2:nq1) = n%balpha(1:nq-2)*venergy(1:nq-2)/units%c**2 ![g/s]
      bradcop(2:nq1) = n%balpha(2:nq1)*o%brc(2:nq1)
     &  * bhrad(2:nq1)*o%bda(2:nq1)                               ![g/s]
      bvye(2:nq1) = n%balpha(1:nq-2)*vye(1:nq-2)                  ![g/s]

!.....advection on the a-grid...........................................
      aurel(1:nq) = -n%ax(1:nq)/dt                                ![g/s]
      btmp(1:nq) = 1./n%bD(1:nq)                               ![cm^3/g]
      call atvdflux(nq,n%ada,n%bda,aurel,btmp,Dflux)           ![cm^3/s]
      call atvdflux(nq,n%ada,n%bda,aurel,n%bei,aEflux)            ![g/s]
      call atvdflux(nq,n%ada,n%bda,aurel,n%bin%ye,yeflux)         ![g/s]
      do inu=1,2                                                  ![g/s]
        call atvdflux(nq,n%ada,n%bda,aurel,n%bynu(1,inu),ynuflux(1,inu))
        call atvdflux(nq,n%ada,n%bda,aurel,n%bznu(1,inu),znuflux(1,inu))
      enddo                                                     ![?*g/s]

!.....advection on the b-grid...........................................
      burel(1) = aurel(1)
      burel(2:nq1) = 0.5*(aurel(2:nq1) + aurel(1:nq-2))           ![g/s]
      burel(nq) = aurel(nq)
      call btvdflux(nq,n%ada,n%bda,burel,n%aS,Sflux)              ![g/s]
      call btvdflux(nq,n%ada,n%bda,burel,n%aek,bEflux)            ![g/s]
      btmp(1) = n%aep(1)
      btmp(2:nq1) = 0.5*(n%aep(2:nq1) + n%aep(1:nq-2))               ![]
      btmp(nq) = n%aep(nq)
      bEflux(1:nq) = bEflux(1:nq) + btmp(1:nq)*burel(1:nq)        ![g/s]

!.....set scale for artificial viscosity................................
      interior = .false.
      interior(2:nq1) = .true.
      eta(2:nq) = 2.*n%ar(2:nq)*input_allocation%visc              ![cm]

!.....artificial viscosity..............................................
      where (interior)
        divu = bdw/bdV                                           ![1/cm]
        comp = bdu/bdr - divu/3.                                 ![1/cm]
      else where
        divu = 0.
        comp = 0.
      end where
      where (divu.lt.0.)
        bQ = eta**2 * n%bin%n*divu*comp*units%c**2           ![g/cm/s^2]
      else where
        bQ = 0.
      end where
      baVQ = n%balpha*bV*bQ                                ![g*cm^2/s^2]
      aQ(1) = bQ(1)
      aQ(2:nq1) = 0.5*(bQ(3:nq) + bQ(2:nq1))                 ![g/cm/s^2]
      aQ(nq) = bQ(nq)
       
!-----calculate equations-----------------------------------------------

!.....center boundaries.................................................
      f(jf(1,9)) = n%ax(1)                                          ![g]
      f(jf(1,10)) = n%ar(1) - o%ar(1)                              ![cm]
      f(jf(1,11)) = n%au(1) - o%au(1)                                ![]
      f(jf(1,12)) = n%am(1) - o%am(1)                               ![g]

!.....center ghost zones................................................
      f(jf(2,1)) = n%bin%n(2) - n%bin%n(1)                     ![g/cm^3]
      f(jf(2,2)) = n%bin%t(2) - n%bin%t(1)                        ![MeV]
      f(jf(2,3)) = n%balpha(2) - n%balpha(1)                         ![]
      f(jf(2,4)) = n%bin%ye(2) - n%bin%ye(1)                         ![]
      f(jf(2,5:6)) = n%bynu(2,1:2) - n%bynu(1,1:2)                   ![]
      f(jf(2,7:8)) = n%bznu(2,1:2) - n%bznu(1,1:2)                   ![]

!.....continuity equation...............................................
      advect(2:nq) = Dflux(2:nq) - Dflux(1:nq1)                ![cm^3/s]
      atmp(1:nq) = -aalpha(1:nq)*n%aw(1:nq)*units%c            ![cm^3/s]
      if (input_idsa%hydro.eq.1) then
        grad(2:nq) = atmp(2:nq) - atmp(1:nq1)                  ![cm^3/s]
      else
        grad(2:nq) = 0.
      endif
      dDdt(2:nq) = (bdV(2:nq)-obdV(2:nq))/dt + advect(2:nq)    ![cm^3/s]
      f(jf(3:nz1,1)) = dDdt(2:nq1) + grad(2:nq1)               ![cm^3/s]

!.....total energy equation.............................................
      atmp(2:nq1) = bEflux(3:nq) - bEflux(2:nq1)                  ![g/s]
      btmp(2:nq1) = aEflux(2:nq1) - aEflux(1:nq-2)                ![g/s]
      where (interior)
        dEdt = (n%bei*n%bda - o%bei*o%bda)/dt + btmp
     1  + ( (n%aek+n%aep)*n%ada - (o%aek+o%aep)*o%ada )/dt + atmp ![g/s]
        btmp = n%balpha*(0.5*(bp + o%bp) + bQ)               ![g/cm/s^2]
        source = -o%ag*bvenergy - o%au*avstress - bradcop         ![g/s]
      end where
      btmp(nq) = n%balpha(nq)*(0.5*(bp(nq)+opsurf)+bQ(nq))   ![g/cm/s^2]
      if (input_idsa%hydro.eq.1) then
        grad(2:nq1) = (n%aw(2:nq1)*btmp(3:nq)
     &    - n%aw(1:nq-2)*btmp(2:nq1))/units%c                     ![g/s]
      else                                                  
        grad(2:nq1) = 0.
      endif
      f(jf(3:nz1,2)) = dEdt(2:nq1) + grad(2:nq1) + source(2:nq1)  ![g/s]

!.....momentum equation.................................................
      advect(1:nq1) = Sflux(2:nq) - Sflux(1:nq1)                  ![g/s]
      atmp(1:nq1) = n%aV(1:nq1)
     1  * (n%bg(2:nq)*n%balpha(2:nq)*bp(2:nq)
     2  - n%bg(1:nq1)*n%balpha(1:nq1)*bp(1:nq1))
     3  + (n%bg(2:nq)*baVQ(2:nq)-n%bg(1:nq1)*baVQ(1:nq1))  ![g*cm^2/s^2]
      apq(2:nq1) = 0.5*(bp(3:nq) + bp(2:nq1)) + aQ(2:nq1)    ![g/cm/s^2]
      an(2:nq1) = 0.5*(n%bin%n(3:nq) + n%bin%n(2:nq1))         ![g/cm^3]
      ae(2:nq1) = 0.5*(n%bout%e(3:nq) + n%bout%e(2:nq1))             ![]
      btmp(2:nq1) = bp(2:nq1)*bjrad(2:nq1)
     &  + bkrad(2:nq1)*(1.+n%bout%e(2:nq1))                  ![g/cm/s^2]
      btmp(nq) = btmp(nq1)
      arad(2:nq1) = 0.5*(btmp(3:nq) + btmp(2:nq1))           ![g/cm/s^2]
      if (input_idsa%hydro.eq.1) then
        if (input_action%rel) then
          where (interior)
            dSdt = (n%aS*n%ada - o%aS*o%ada)/dt + advect          ![g/s]
            grad = 3.*atmp/(n%ar*units%c)                         ![g/s]
            source =aalpha*( (1.+ae)*n%amr*(units%c**2+6.*n%aV*apq/n%am)
     1    + 4.*units%pi*n%ar**2 * units%G/units%c**2*arad
     2    + (n%au**2*(2.*apq-3.*aQ)-n%amr*apq)/an )/(o%ar*units%c)*n%ada
     &    - o%ag*avstress - o%au*bvenergy                         ![g/s]
          end where
        else
          where (interior)
            dSdt = (n%aS*n%ada - o%aS*o%ada)/dt + advect          ![g/s]
            grad = 3.*atmp/(n%ar*units%c)                         ![g/s]
            source = n%amr*units%c/o%ar*n%ada - avstress          ![g/s]
          end where
        endif
      else
        where (interior)
          dSdt = (n%aS*n%ada - o%aS*o%ada)/dt + advect            ![g/s]
          grad = 0.
          source = - o%ag*avstress - o%au*bvenergy                ![g/s]
        end where
      endif
      f(jf(3:nz1,3)) = dSdt(2:nq1) + grad(2:nq1) + source(2:nq1)  ![g/s]

!.....electron fraction equation........................................
      advect(2:nq) = yeflux(2:nq) - yeflux(1:nq1)                 ![g/s]
      where (interior)
        dyedt = (n%bin%ye*n%bda - o%bin%ye*o%bda)/dt + advect     ![g/s]
        source = -bvye                                            ![g/s]
      end where
      f(jf(3:nz1,4)) = dyedt(2:nq1) + source(2:nq1)               ![g/s]

!.....neutrino fraction equation........................................
      do inu=1,2
        advect(2:nq) = ynuflux(2:nq,inu) - ynuflux(1:nq1,inu)     ![g/s]
        do i=2,nq1
          dnudt(i) = (n%bynu(i,inu)*n%bda(i)                      ![g/s]
     &       - o%bynu(i,inu)*o%bda(i))/dt + advect(i)             ![g/s]
          source(i) = -n%balpha(i)*vynu(inu,i-1)                  ![g/s]
        enddo
        f(jf(3:nz1,4+inu)) = dnudt(2:nq1) + source(2:nq1)         ![g/s]
        advect(2:nq) = znuflux(2:nq,inu) - znuflux(1:nq1,inu)   ![?*g/s]
        do i=2,nq1
          dnudt(i) = (n%bznu(i,inu)*n%bda(i)
     &       - o%bznu(i,inu)*o%bda(i))/dt + advect(i)           ![?*g/s]
          source(i) = -n%balpha(i)*vznu(inu,i-1)                ![?*g/s]
        enddo
        f(jf(3:nz1,6+inu)) = dnudt(2:nq1) + source(2:nq1)       ![?*g/s]
      enddo

!.....volume integration................................................
      f(jf(2:nz1,10)) = n%bda(2:nq) - n%bD(2:nq)*bdV(2:nq)          ![g]

!.....gravitational mass integration....................................
      btmp(2:nq) = n%am(2:nq) - n%am(1:nq1)                         ![g]
      if (input_action%rel) then
        f(jf(2:nz1,11)) = btmp(2:nq) - (n%bg(2:nq)+n%bei(2:nq)      ![g]
     1    + n%ag(2:nq)*bjrad(2:nq) + n%au(2:nq)*bhrad(2:nq))*n%bda(2:nq)
      else
        f(jf(2:nz1,11)) = btmp(2:nq) - n%bda(2:nq)                  ![g]
      endif

!.....lapse function....................................................
      atmp(2:nq1) = an(2:nq1)*(1. + ae(2:nq1))*units%c**2    ![g/cm/s^2]
      if (input_action%rel) then
        f(jf(2:nz-2,12)) =
     1    (n%balpha(3:nq)*bp(3:nq) - n%balpha(2:nq1)*bp(2:nq1))
     2    + atmp(2:nq1)*(n%balpha(3:nq) - n%balpha(2:nq1))
     3    + (baVQ(3:nq) - baVQ(2:nq1))/n%aV(2:nq1)           ![g/cm/s^2]
     4    - avstress(2:nq1)*units%c/(4.*units%pi*n%ar(2:nq1)**2)
      else
        f(jf(2:nz-2,12)) = n%balpha(3:nq) - n%balpha(2:nq1)          ![]
      endif

!.....surface ghost zones...............................................
      f(jf(nz,1)) = n%bin%n(nq) - o%bin%n(nq)                  ![g/cm^3]
      f(jf(nz,2)) = n%bin%t(nq) - o%bin%t(nq)                     ![MeV]
      f(jf(nz,3)) = aalpha(nq)*n%aw(nq) - aalpha(nq1)*n%aw(nq1)  ![cm^2]
      f(jf(nz,4)) = n%bin%ye(nq) - o%bin%ye(nq)                      ![]
      f(jf(nz,5:6)) = n%bynu(nq,1:2) - o%bynu(nq,1:2)                ![]
      f(jf(nz,7:8)) = n%bznu(nq,1:2) - o%bznu(nq,1:2)               ![?]

!.....surface boundaries................................................
      f(jf(nz1,9)) = n%ax(nq1)                                      ![g]
      if (input_action%rel) then
        f(jf(nz1,12)) = n%balpha(nq) - (1.-2.*n%amr(nq))/n%ag(nq)    ![]
      else
        f(jf(nz1,12)) = n%balpha(nq) - 1.                            ![]
      endif

      end subroutine equations_step

!=======================================================================
