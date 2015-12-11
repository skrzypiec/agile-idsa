!=======================================================================
!
!     IDSA: calculates neutrino approximations in spherical symmetry
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine nuprox(n,dt,dV,d,s,ye,ynu,znu,iter,iray,
     &  sout,edot,yout,yedot,vdot,ynudot,znudot,dLmax)

      use egroup_module
      use input_file_module
      use mu_tau_leakage_module
      use species_module
      use spectrum_module
      use units_module

      implicit none

      integer, intent(in) :: n,iray
      real, intent(in) :: dt
      real, dimension(n), intent(in) :: dV,d,s,ye
      real, dimension(nut,n), intent(in) :: ynu,znu
      integer, intent(inout) :: iter
      real, dimension(n), intent(out) :: sout,edot,yout,yedot,vdot
      real, dimension(nut), intent(out) :: dLmax
      real, dimension(nut,n), intent(out) :: ynudot,znudot

!-----------------------------------------------------------------------
!
!     Implements the isotropic diffusion source approximation for
!     neutrino transport in a spherically symmetric supernova model
!     along Liebendoerfer, Whitehouse & Fischer, ApJ 698 (2009).
!
!     Input:
!     n ... number of shells
!     dt ... time step [s]
!     dV ... zone volume [cm^3], dV(1) is a sphere around the origin
!     d  ... baryon mass density [g/cm^3]
!     s  ... specific entropy [kB/baryon]
!     t  ... matter temperature [MeV]
!     ye ... electron fraction []
!     ynu ... neutrino fraction []
!     znu ... neutrino energy [MeV/baryon]
!     iter ... maximum number of iterations
!     iray ... 1 for a 1-D call, 2 for a 3-D call
!
!     Output:
!     iter ... number of iterations taken
!     sout ... specific entropy after time step [kB/baryon]
!     edot ... specific internal energy change rate [erg/g/s]
!     yout ... electron fraction after time step []
!     yedot ... electron fraction change rate [1/s]
!     vdot ... outer zone edge neutrino stress [cm/s^2]
!     ynudot ... neutrino fraction change rate [1/s]
!     znudot ... neutrino energy change rate [MeV/baryon/s]
!     dLmax ... relative change in number luminosity at infinity
!
!     Static in spectral_module:
!     src ... spectral neutrino source, sigma-chi*fs [particles/baryon/s]
!
!     Note: zone n is not updated
!
!-----------------------------------------------------------------------

      logical, parameter :: implicit_diffusion=.true.
      real, parameter :: cvr=2./3.
      real, parameter :: nuproplength=1.5e+06
      real, parameter :: fac=4.*units%pi/(units%h*units%c)**3

      integer, dimension(:), allocatable :: errflag
      integer :: maxiter,status,i,it,ie
      real :: tmp,ff,xcoff,ycoff,zcoff,vol
      real :: qloc,qdiff,qeff,epsilon
      real, dimension(ne) :: ps
      real, dimension(:), allocatable :: t,ar,br,apir2,bpir2,dm,pnu
      real, dimension(nut) :: ynutmp,znutmp
      real, dimension(ne,nut) :: feq,etai,alpha,ft,rnu
      real, dimension(:,:,:), allocatable :: mfpe,f,dfdt,taue,fs
      real, dimension(:,:,:), allocatable :: flux,sl,amfpe,srccp
      real, dimension(:,:), allocatable :: mfpm,amfpm,taum
      real, dimension(:,:), allocatable :: firstL,Lold,Lnew

!.....allocate variables................................................
      allocate(errflag(0:n+1))
      allocate(t(n),ar(n),br(n),apir2(n),bpir2(n),dm(n),pnu(n))
      allocate(mfpe(ne,nut,n),mfpm(ne,n),f(ne,nut,n),dfdt(ne,nut,n),
     &  taue(ne,nut,n),taum(ne,n),fs(ne,nut,n),srccp(ne,nut,n))
      allocate(flux(ne,nut,n),sl(ne,nut,n),amfpe(ne,nut,n),amfpm(ne,n))
      allocate(firstL(ne,nut),Lold(ne,nut),Lnew(ne,nut))

!.....distances (a=zone edge value, b=zone center value)................
      dm(1) = d(1)*dV(1)
      vol = dV(1)
      br(1) = 0.                                        !origin
      ar(1) = (3.*vol/(4.*units%pi))**(1./3.)           !radius of first zone
      apir2(1) = 4.*units%pi*ar(1)*ar(1)                !cm^2
      bpir2(1) = 0.                                     !cm^2 (not used)
      do i=2,n
        dm(i) = d(i)*dV(i)                              !g
        vol = vol + 0.5*dV(i)
        br(i) = (3.*vol/(4.*units%pi))**(1./3.)         !cm
        vol = vol + 0.5*dV(i)
        ar(i) = (3.*vol/(4.*units%pi))**(1./3.)         !cm
        apir2(i) = 4.*units%pi*ar(i)*ar(i)              !cm^2
        bpir2(i) = 4.*units%pi/3.
     &    * (ar(i)**2+ar(i)*ar(i-1)+ar(i-1)**2)         !cm^2
      enddo

!.....phase space.......................................................
      if (nzone.lt.n) then
        write(6,*) 'Error: nzone smaller than n = ',n
        stop
      endif
      if (egroup_empty) call egroup
      ps = fac*E**2*dE*units%mb	              !particles/cm^3 * g/baryon

!.....iterate for convergence...........................................
      maxiter = iter
      epsilon = input_step%eps
      do iter=1,maxiter

!.....copy sources for luminosity output by wrnuprox()..................
        srccp = src(:,:,1:n,iray)

!-----solving the Poisson equation--------------------------------------

!.....calculate particle flux on zone edges.............................
        do it=1,nut
          do ie=1,ne
            tmp = 0.
            do i=1,n
              tmp = tmp + src(ie,it,i,iray)*dm(i)/units%mb
              flux(ie,it,i) = tmp/apir2(i)    !particles/s/cm^2
            enddo
          enddo
        enddo

!-----finding the geometry of the neutrino scattering spheres-----------

!.....determine distribution function and mean free path................
        do i=1,n
          call meanfreepath(i,d(i),s(i),ye(i),ynu(1,i),znu(1,i),iray,
     &      t(i),f(1,1,i),feq,mfpe(1,1,i),mfpm(1,i),status)
          errflag(i) = status
        enddo
        errflag(0) = errflag(2)     !mirror at origin for update (0,1,2)
        errflag(n+1) = 16           !disable update (n-1,n,n+1)

!.....calculate optical depth for electron flavour neutrinos............
        amfpe(:,:,n) = 1.d+99    !mean free path on zone edge [cm]
        taue(:,:,n) = 0.         !optical depth on zone edge []
        do it=1,nut
          do ie=1,ne
            do i=n,2,-1
              if (errflag(i).eq.0) then
                amfpe(ie,it,i-1) = 2./(1./mfpe(ie,it,i)
     &            + 1./mfpe(ie,it,i-1))   
                taue(ie,it,i-1) = taue(ie,it,i)
     &            + (ar(i)-ar(i-1))/mfpe(ie,it,i)
              else
                amfpe(ie,it,i-1) = amfpe(ie,it,i)
                taue(ie,it,i-1) = taue(ie,it,i)
              endif
            enddo

!.....find neutrinospheres.............................................
            if (taue(ie,it,1).le.cvr) then
              rnu(ie,it) = 0.           !whole domain transparent
            else
              i = n
              do while (taue(ie,it,i).le.cvr)
                i = i-1
              enddo
              rnu(ie,it) = ar(i) + (ar(i+1)-ar(i))
     &          * (cvr-taue(ie,it,i))/(taue(ie,it,i+1)-taue(ie,it,i))   !cm
            endif
          enddo
        enddo

!.....calculate optical depth for mu/tau flavour neutrinos..............
        amfpm(:,n) = 1.d+99      !mean free path on zone edge [cm]
        taum(:,n) = 0.           !optical depth on zone edge []
        do ie=1,ne
          do i=n,2,-1
            if (errflag(i).eq.0) then
              amfpm(ie,i-1) = 2./(1./mfpm(ie,i) + 1./mfpm(ie,i-1))
              taum(ie,i-1) = taum(ie,i)
     &          + (ar(i)-ar(i-1))/mfpm(ie,i)
            else
              amfpm(ie,i-1) = amfpm(ie,i)
              taum(ie,i-1) = taum(ie,i)
            endif
          enddo
        enddo

!-----implicit diffusion update-----------------------------------------

!.....proceed with the remaining valid zones............................
        do i=1,n

!.....estimate the streaming particle density...........................
          do it=1,nut
            do ie=1,ne
              if (i.eq.1) then
                fs(ie,it,i) = 0.     !no streaming particles at center
              else
                if (abs(rnu(ie,it)).ge.br(i)) then
                  ff = 0.5
                else
                  ff = 0.5*(1. + sqrt(1.-(rnu(ie,it)/br(i))**2))
                endif
                fs(ie,it,i) = units%mb*max(flux(ie,it,i-1),0.)
     &            * apir2(i-1)/(d(i)*bpir2(i)*units%c*ff)
                                              !particles/baryon
              endif
            enddo
          enddo

!.....check validity of zone and its neighbors..........................
          if (errflag(i).eq.0) then
            if (   ((errflag(i-1).eq.0).or.(errflag(i-1).gt.20))
     &        .and.((errflag(i+1).eq.0).or.(errflag(i+1).gt.20)) ) then
              do it=1,nut
                do ie=1,ne

!.....calculate the diffusion term......................................
                  tmp = units%c/(3.*dV(i))                     !1/s/cm^2
                  xcoff = tmp*apir2(i)*amfpe(ie,it,i)/(br(i+1)-br(i))
                  if (i.eq.1) then
                    etai(ie,it) = 0.                                !1/s
                    alpha(ie,it) = xcoff*(f(ie,it,i)-f(ie,it,i+1))  !1/s
                  else
                    zcoff = tmp*apir2(i-1)*amfpe(ie,it,i-1)
     &                    / (br(i)-br(i-1))
                    ycoff = 0.5*(xcoff+zcoff)
                    if (implicit_diffusion) then
                      etai(ie,it) = zcoff                           !1/s
                    else
                      etai(ie,it) = 0.
                    endif
                    alpha(ie,it) = 2.*ycoff*f(ie,it,i)
     &                - xcoff*f(ie,it,i+1) - zcoff*f(ie,it,i-1)     !1/s
                  endif
                enddo
                ft(:,it) = ps*f(:,it,i)/d(i)           !particles/baryon
                alpha(:,it) = ps*alpha(:,it)/d(i)    !particles/baryon/s
              enddo

!.....calculate new sources.............................................
              call uyupdate(dt,d(i),s(i),ye(i),nuproplength,
     &          ft,fs(1,1,i),etai,alpha,iray,
     &          sout(i),yout(i),dfdt(1,1,i),sl(1,1,i),status)
              if (status.ne.0) errflag(i) = status
            else
              errflag(i) = 15          !zone has invalid neighbour zones
            endif
          endif

!.....set dummy updates if diffusion step failed........................
          if (errflag(i).ne.0) then
            sout(i) = s(i)
            yout(i) = ye(i)
            dfdt(:,:,i) = 0.
            sl(:,:,i) = 0.
          endif

!.....calculate change rates............................................
          src(:,:,i,iray) = sl(:,:,i) - dfdt(:,:,i)
          yedot(i) = -sum(sl(:,len,i)-sl(:,lea,i))   !1/s
          edot(i) = -sum((sl(:,len,i)+sl(:,lea,i))*E)
     &      * units%MeV/units%mb	             !erg/g/s
          do it=1,nut
            ynudot(it,i) = sum(dfdt(:,it,i))         !particles/baryon/s
            znudot(it,i) = sum(dfdt(:,it,i)*E)       !MeV/baryon/s
          enddo

!.....leakage scheme for the cooling by emission of mu/tau neutrinos....
          if ((input_idsa%mutau.eq.1).and.(errflag(i).eq.0)) then
            call leakage(d(i),t(i),ye(i),mfpm(1,i),taum(1,i),
     &        qloc,qdiff,qeff)
!           write(6,11) iray,i,d(i),t(i),ye(i),edot(i),-qeff/d(i)
11          format(i6,i6,5g12.4)
            edot(i) = edot(i) - qeff/d(i)
          endif

!.....update spectrum...................................................
          call spectrum_update(i,dt,d(i),ynu(:,i),znu(:,i),dfdt(:,:,i),
     & iray)

!.....update distribution function......................................
          ynutmp = ynu(:,i) + ynudot(:,i)*dt
          znutmp = znu(:,i) + znudot(:,i)*dt
          pnu(i) = d(i)/units%mb*sum(znutmp)*units%MeV/3.
          if (implicit_diffusion) then
            call meanfreepath(i,d(i),sout(i),yout(i),ynutmp,znutmp,iray,
     &        t(i),f(1,1,i),feq,mfpe(1,1,i),mfpm(1,i),status)
          endif
          if (status.ne.0) errflag(i) = status
        enddo

!.....calculate luminosity changes......................................
        Lold = 0.
        Lnew = 0.
        do it=1,nut
          do i=1,n
            Lold(:,it) = Lold(:,it)+srccp(:,it,i)*dm(i)/units%mb
            Lnew(:,it) = Lnew(:,it)+src(:,it,i,iray)*dm(i)/units%mb
          enddo
        enddo
        if (iter.eq.1) firstL = Lold
        dLmax = 0.
        do it=1,nut
          tmp = 0.
          do ie=1,ne
            tmp = max(tmp,abs(Lnew(ie,it)))
            dLmax(it) = max(dLmax(it),abs(Lnew(ie,it)-Lold(ie,it)))
          enddo
          dLmax(it) = dLmax(it)/(tmp+1.d-38)
        enddo
        if ((dLmax(len).lt.epsilon).and.(dLmax(lea).lt.epsilon)) exit
      enddo
      dLmax = 0.
      do it=1,nut
        tmp = 0.
        do ie=1,ne
          tmp = max(tmp,abs(Lnew(ie,it)))
          dLmax(it) = max(dLmax(it),abs(Lnew(ie,it)-firstL(ie,it)))
        enddo
        dLmax(it) = dLmax(it)/(tmp+1.d-38)
      enddo

!.....calculate neutrino stress.........................................
      do i=1,n-1
        vdot(i) = -2.*apir2(i)*(pnu(i+1)-pnu(i))/(dm(i+1)+dm(i)) ![cm/s^2]
      enddo
      vdot(n) = 0.

!.....write result......................................................
      call wrnuprox(n,dV,d,s,ye,ynu,znu,fs,srccp,sl,mfpe,iray,
     &  errflag)

      deallocate(errflag)
      deallocate(t,ar,br,apir2,bpir2,dm,pnu,firstL,Lold,Lnew)
      deallocate(mfpe,mfpm,f,dfdt,taue,taum,fs,flux,sl,amfpe,srccp)

      end subroutine nuprox

!=======================================================================
