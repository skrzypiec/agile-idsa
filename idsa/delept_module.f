!***********************************************************************
!
!     delept_module, cheap handling of deleptonization during collapse
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module delept_module

      use units_module

      implicit none

      real, parameter :: delept_dtrap = 2.e+12  !trapping dens. [g/cm^3]
      real, parameter :: delept_Enu = 10.       !escape energy [MeV]
      real, parameter :: delept_DeltaQ = units%Q!n-p mass diff. [MeV]

!-----------------------------------------------------------------------
!
!     This module estimates the Ye and entropy changes during the
!     collapse phase as described in Liebendoerfer, ApJ 633:1042 (2005).
!
!-----------------------------------------------------------------------

      contains

!=======================================================================
!
!     delept_ysv, interface for parmeterised deleptonisation
!
!=======================================================================

        subroutine delept_ysv(n,dt,dV,d,s,ye,edot,ydot,vdot)

        use eos_module

        integer, intent(in) :: n
        real, intent(in) :: dt
        real, dimension(n), intent(in) :: dV,d,s,ye
        real, dimension(n), intent(out) :: edot,ydot,vdot

!-----------------------------------------------------------------------
!
!     Input:
!     n ... number of shells
!     dt ... time step [s]
!     dV ... zone volume [cm^3], dV(1) is a sphere around the origin
!     d  ... baryon mass density [g/cm^3]
!     s  ... specific entropy [kB/baryon]
!     ye ... electron fraction []
!
!     Output:
!     edot ... specific internal energy change rate [erg/g/s]
!     yedot ... electron fraction change rate [1/s]
!     vdot ... outer zone edge neutrino stress [cm/s^2]
!
!-----------------------------------------------------------------------

        integer :: i,status
        real :: sdot,cs,e0,y1,s1
        type(eos_state) :: state

!.....parameterisation scheme for Ye and entropy........................
        do i=1,n
          call delept_ys(dt,d(i),s(i),ye(i),ydot(i),sdot)

!.....convert entropy change to energy change...........................
          call eos_sinterp(d(i),ye(i),s(i),cs,state,status)
          if (status.ne.0) then
            write(6,*) 'Warning: status = ',status,' in agile_delept'
            edot(i) = 0.
            goto 100
          endif
          e0 = state%e                                          ![erg/g]
          y1 = ye(i)+ydot(i)*dt
          s1 = s(i)+sdot*dt
          call eos_sinterp(d(i),y1,s1,cs,state,status)
          if (status.ne.0) then
            write(6,*) 'Warning: status = ',status,' in agile_delept'
            edot(i) = 0.
            goto 100
          endif
          edot(i) = (state%e-e0)/dt                           ![erg/g/s]
100       continue
        enddo

!.....parameterisation scheme for neutrino stress.......................
        call delept_v(n,dt,dV,d,s,ye,ydot,vdot)

        end subroutine delept_ysv

!=======================================================================
!
!     delept_ys, estimates the Ye and entropy changes
!
!=======================================================================

        subroutine delept_ys(dt,d,s,ye,ydot,sdot)

        use species_module

        real, intent(in) :: dt,d,ye,s
        real, intent(out) :: ydot,sdot

!-----------------------------------------------------------------------
!
!     Input:
!     dt ... time step [s]
!     d   ... density [g/cm^3]
!     s   ... entropy [kB/baryon]
!     ye  ... electron fraction []
!     parameter ::
!    &  logdi=(/7.47712125471966,13.30102999566398/) ...density [g/cm^3]
!     parameter :: yi=(/ 0.5, 0.285 /)          ... electron fraction []
!     real, parameter :: yc=0.035
!
!     Output:
!     ydot ... electron fraction change rate [1/s]
!     sdot  ... entropy change [kB/baryon/s]
!
!-----------------------------------------------------------------------

        real, dimension(2), parameter ::
     &    logdi=(/7.30102999566398,13.30102999566398/)
        real, dimension(2), parameter :: yi=(/ 0.5, 0.285 /)
        real, parameter :: yc=0.035
        integer :: imaxm,j,jl,jm,ju,status
        real :: logd,tmp,x,ax,lambda,ytmp,u,munu,t
        real, dimension(nuc+1) :: y,mu

!.....calculate electron fraction change................................
        logd = log10(d)

!.....reference Ye given as fitting formula.............................
        tmp = (2.*logd-logdi(2)-logdi(1))/(logdi(2)-logdi(1))
        x = max(-1.,min(1.,tmp))
        ax = abs(x)
        ytmp = 0.5*(yi(2)+yi(1))+0.5*x*(yi(2)-yi(1))
     &    + yc*( 1.-ax+4.*ax*(ax-0.5)*(ax-1.) )

!.....call equation of state............................................
        status = 0
        y(je) = ye
        call eosinterf(d,s,y,mu,t,u,status)

!.....calculate Ye and entropy changes..................................
        ydot = min(0.d0,ytmp-ye)/dt
11      format(4g12.4)
        sdot = 0.
        if (status.eq.0) then
          munu = mu(je) - mu(jhat) - delept_DeltaQ ![MeV]
          if ((munu.gt.delept_Enu).and.(d.lt.delept_dtrap))
     &      sdot = -ydot*(munu-delept_Enu)/t
        endif

        end subroutine delept_ys

!=======================================================================
!
!     delept_v, determines the neutrino stress based on the luminosity
!
!=======================================================================

        subroutine delept_v(n,dt,dV,d,s,ye,ydot,vdot)

        use species_module

        integer, intent(in) :: n
        real, intent(in) :: dt
        real, dimension(n), intent(in) :: dV,d,ye,s,ydot
        real, dimension(n), intent(out) :: vdot

!-----------------------------------------------------------------------
!
!     Input:
!     n ... number of shells
!     dt ... time step [s]
!     dV ... zone volume [cm^3], dV(1) is a sphere around the origin
!     d  ... baryon mass density [g/cm^3]
!     s  ... specific entropy [kB/baryon]
!     ye ... electron fraction []
!     ydot ... electron fraction change rate [1/s]
!
!     Output:
!     vdot  ... outer zone edge neutrino stress [cm/s^2]
!
!-----------------------------------------------------------------------

        logical :: shocked
        integer :: status,i,i0
        real, parameter :: pi=units%pi
        real :: lsum,u,munu,eta,fak,vol
        real, dimension(nuc+1) :: y,mu
        real, dimension(:), allocatable :: lum,t,pnu,ar,br,apir2,dm
        real, dimension(:), allocatable :: F2,F3,F5

!.....allocate variables................................................
        allocate(lum(n),t(n),pnu(n),ar(n),br(n),apir2(n),dm(n),
     &    F2(0:n),F3(0:n),F5(0:n))

!.....distances (a=zone edge value, b=zone center value)................
        dm(1) = d(1)*dV(1)
        vol = dV(1)
        br(1) = 0.                                        !origin
        ar(1) = (3.*vol/(4.*units%pi))**(1./3.)    !radius of first zone
        apir2(1) = 4.*units%pi*ar(1)*ar(1)                !cm^2
        do i=2,n
          dm(i) = d(i)*dV(i)                              !g
          vol = vol + 0.5*dV(i)
          br(i) = (3.*vol/(4.*units%pi))**(1./3.)         !cm
          vol = vol + 0.5*dV(i)
          ar(i) = (3.*vol/(4.*units%pi))**(1./3.)         !cm
          apir2(i) = 4.*units%pi*ar(i)*ar(i)              !cm^2
        enddo

!.....Fermi integrals at eta=0............................................
        F2(0) = 3./2.*1.202
        F3(0) = 7.*pi**4./120.
        F5(0) = 31.*pi**6./252.

        lsum = 0.
        do i=1,n

!.....number luminosity.................................................
          lsum = lsum - dm(i)/units%mb*ydot(i) ![particles/s]
          lum(i) = lsum/apir2(i) ![particles/cm^2/s]

!.....call equation of state............................................
          status = 0
          y(je) = ye(i)
          call eosinterf(d(i),s(i),y,mu,t(i),u,status)

!.....neutrino pressure.................................................
          if (status.eq.0) then
            munu = mu(je) - mu(jhat) - delept_DeltaQ ![MeV]
            eta = munu/t(i)                                  ![ ]
            if (eta.gt.0.) then
              F2(i) = (eta**3. + pi**2.*eta)/3.  
     &              + 3./2.*1.202*exp(-0.825*eta)
              F3(i) = (eta**4. + 2.*pi**2.*eta**2. + 7.*pi**4./15.)/4.  
     &              - 7.*pi**4./120.*exp(-eta)
              F5(i) = (eta**6. + 5.*pi**2*eta**4. + 7.*pi**4.*eta**2.  
     &              + 31.*pi**6./21.)/6. - 31.*pi**6./252.*exp(-eta)
            else
              F2(i) = 3./2.*1.202*exp(eta)
              F3(i) = 7.*pi**4./120.*exp(eta)
              F5(i) = 31.*pi**6./252.*exp(eta)
            endif
            pnu(i) = 4.*pi/3.*units%MeV/(units%h*units%c)**3.  
     &        * t(i)**4. * F3(i) ![erg/cm^3]
          else
            F2(i) = F2(0)
            F3(i) = F3(0)
            F5(i) = F5(0)
            if (i.eq.1) then
              pnu(i) = 0.
            else
              pnu(i) = pnu(i-1)
            endif
          endif
        enddo

!.....neutrino stress...................................................
        vdot = 0.
        i = 1
        do while (d(i).gt.delept_dtrap)
          vdot(i) = -2.*apir2(i)*(pnu(i+1)-pnu(i))/(dm(i+1) + dm(i)) ![cm/s^2]
          i = i+1
        enddo
        i0 = i-1
        if (i0.ge.1) then
          if (s(i0).lt.3) then
            fak = vdot(i0)  
     &        / (t(i0+1)**3.*F5(i0+1)/F2(i0+1)*lum(i0))
            do i=i0+1,n-1
              vdot(i) = fak*max(t(i+1)**3.*F5(i+1)/F2(i+1),  
     &          delept_Enu**3.*F2(0)**2.*F5(0)/F3(0)**3.)  
     &          * lum(i) ![cm/s^2]
            enddo
          endif
        endif

!.....deallocate variables..............................................
        deallocate(lum,t,pnu,ar,br,apir2,dm,F2,F3,F5)

        end subroutine delept_v

!=======================================================================

      end module delept_module

!***********************************************************************
