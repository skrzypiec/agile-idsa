!=======================================================================
!
!     IDSA: wrnuprox
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine wrnuprox(n,dV,d,s,ye,ynu,znu,fs,srccp,sl,lambda,iray,
     &  errflag)

      use agile_module, only: istep
      use egroup_module
      use input_file_module
      use species_module
      use spectrum_module
      use units_module

      implicit none

      integer, intent(in) :: n,iray
      integer, dimension(0:n+1) :: errflag
      real, dimension(n), intent(in) :: dV,d,s,ye
      real, dimension(nut,n), intent(in) :: ynu,znu
      real, dimension(ne,nut,n), intent(in) :: fs,srccp,sl,lambda

!-----------------------------------------------------------------------
!
!     Input:
!     n ... number of shells
!     errflag ... error flag for each zone
!     dV ... zone volume [cm^3], dV(1) is a sphere around the origin
!     d  ... baryon mass density [g/cm^3]
!     s  ... specific entropy [kB/baryon]
!     ye ... electron fraction []
!     ynu ... neutrino fraction []
!     znu ... neutrino energy [MeV/baryon]
!     fs ... neutrino distribution functions of streaming neutrinos
!          fs = f(E) * 4.*pi/(h*c)**3*E**2*dE*mb/d  ![particles/baryon]
!     srccp ... neutrino source leaving and entering zones [particles/baryon/s]
!     sl ... local neutrino source [particles/baryon/s]
!     iray ... 1 for 1D, 2 for 3D
!
!     List of error flags:
!     
!     eos_module.f:
!      1 ... density out of range
!      2 ... ye out of range
!      3 ... entropy out of range
!      4 ... empty grid point in equation of state table
!     
!     eosinterf.f
!      5 ... abundance lower than zero
!     
!     etafit.f
!     11 ... singular Jacobian
!     12 ... negative eta in degenerate regime
!     13 ... no convergence
!     
!     nuprox.f
!     15 ... zone has invalid neighbor zone
!     16 ... zone is the surface boundary
!     
!     framerate.f
!     21 ... log(entropy) larger than zero
!     
!     uyupdate.f
!     22 ... number of iterations larger than itmax
!     23 ... singular Jacobian
!     
!-----------------------------------------------------------------------

      character(7) :: cycstr
      character(58) :: filename
      character (len=2) :: ftype
      integer, save :: cycsav=-1
      integer :: i,status,inu,ie
      integer, dimension(nut) :: statnut
      real, parameter :: fac=4.*units%pi/(units%h*units%c)**3
      real :: vol,t,u,nb,nd,ed,fermi
      real, dimension(nut) :: beta,eta,tnu,yedot,edot,LN,LE
      real, dimension(nuc+1) :: y,mu
      real, dimension(ne) :: ps,tmp
      real, dimension(ne,nut,n) :: ft,LNk,LEk
      real, dimension(ne,nut) :: localft
      real, dimension(:), allocatable :: ar

!.....write data into files if cycle number has changed.................
      if (iray.eq.1) then
        if (istep.eq.cycsav) return
        cycsav = istep
	ftype = '1d'
      else
        write(6,*) 'Error: invalid iray in wrnuprox.f!'
        stop
      endif

!.....generate cycle string.............................................
      write(cycstr,11) istep,'.d'
11    format(i5,a2)
      do i=1,5
        if(cycstr(i:i).eq.' ') cycstr(i:i)='0'
      enddo
      if (input_path%initial) call input_read_path

!.....open files for output.............................................
      filename = trim(input_path%data)//'nuprox'//ftype//cycstr
      write(6,22) filename
22    format(1x,'Info: writing ',a58)
      open(1,file=filename,status='unknown')
      filename = trim(input_path%data)//'distft'//ftype//cycstr
      open(2,file=filename,status='unknown')
      filename = trim(input_path%data)//'distfs'//ftype//cycstr
      open(3,file=filename,status='unknown')
      filename = trim(input_path%data)//'lambda'//ftype//cycstr
      open(4,file=filename,status='unknown')
      filename = trim(input_path%data)//'LNk'//ftype//cycstr
      open(7,file=filename,status='unknown')
      filename = trim(input_path%data)//'LEk'//ftype//cycstr
      open(8,file=filename,status='unknown')

!.....write headers.....................................................
      write(1,44) 'zone','err','ar','d','T','Ye','mue','muhat',
     &'Ynu','Ynubar','Enu','Enubar','Tnu','Tnubar','etanu','etanubar',
     &'yedot','yedotbar','edot','edotbar','LN','LNbar','LE','LEbar'
      write(1,44) ' ',' ','cm','g/cm^3','MeV',' ','MeV','MeV',
     &' ',' ','MeV','MeV','MeV','MeV',' ',' ',
     &'1/s','1/s','erg/s/g','erg/s/g','1/s','1/s','erg/s','erg/s'
      write(1,55)

      write(2,44) '0','0','E(1)','...','E(kmax)','in MeV'
      write(2,44) '0','0','dE(1)','...','dE(kmax)','in MeV'
      write(2,44) 'zone','err','dYnu(E)','for','trapped',' ',
     &  'neutrinos','in','particles','/baryon/MeV'
      write(2,44) 'zone','err','dYnu(E)','for','trapped',
     &  'anti -','neutrinos','in','particles','/baryon/MeV'
      write(2,55)
      write(2,33) 0,0,E
      write(2,33) 0,0,dE
        
      write(3,44) '0','0','E(1)','...','E(kmax)','in MeV'
      write(3,44) '0','0','dE(1)','...','dE(kmax)','in MeV'
      write(3,44) 'zone','err','dYnu(E)','for','streaming',' ',
     &  'neutrinos','in','particles','/baryon/MeV'
      write(3,44) 'zone','err','dYnu(E)','for','streaming',
     &  'anti -','neutrinos','in','particles','/baryon/MeV'
      write(3,55)
      write(3,33) 0,0,E
      write(3,33) 0,0,dE
      
      write(4,44) '0','0','E(1)','...','E(kmax)','in MeV'
      write(4,44) '0','0','dE(1)','...','dE(kmax)','in MeV'
      write(4,44) 'zone','err','lambda(cm)','for','streaming',' ',
     &  'neutrinos','in','particles','/baryon/MeV'
      write(4,44) 'zone','err','lambda(cm)','for','streaming',
     &  'anti -','neutrinos','in','particles','/baryon/MeV'
      write(4,55)
      write(4,33) 0,0,E
      write(4,33) 0,0,dE

      write(7,44) '0','0','E(1)','...','E(kmax)','in MeV'
      write(7,44) '0','0','dE(1)','...','dE(kmax)','in MeV'
      write(7,44) 'zone','err','LN'
      write(7,55)
      write(7,33) 0,0,E
      write(7,33) 0,0,dE

      write(8,44) '0','0','E(1)','...','E(kmax)','in MeV'
      write(8,44) '0','0','dE(1)','...','dE(kmax)','in MeV'
      write(8,44) 'zone','err','LE'
      write(8,55)
      write(8,33) 0,0,E
      write(8,33) 0,0,dE

44    format(2a6,24a12)
55    format(288('-'))

!.....phase space and luminosity initialisation.........................
      if (egroup_empty) call egroup
      ps = fac*E**2*dE*units%mb !particles/cm^3 * g/baryon
      LN = 0.
      LE = 0.
      LNk = 0.
      LEk = 0.

!.....calculate the outer radius of zones...............................
      allocate(ar(n))
      vol = dV(1)
      ar(1) = (3.*vol/(4.*units%pi))**(1./3.)           !radius of first zone
      do i=2,n
        vol = vol + dV(i)
        ar(i) = (3.*vol/(4.*units%pi))**(1./3.)         !cm
      enddo

!.....loop over zones...................................................
      do i=1,n

!.....call equation of state............................................
        y(je) = ye(i)
        status = 0
        call eosinterf(d(i),s(i),y,mu,t,u,status)
        if (status.ne.0) then
          y = 0.
          mu = 0.
          t = 0.
          u = 0.
        endif

!.....calculate neutrino temperature and degeneracy.....................
        nb = d(i)/units%mb   !baryon density
        do inu=1,nut
          nd = nb*ynu(inu,i)        !neutrino number density
          ed = nb*znu(inu,i)        !neutrino specific energy
          call etafit(nd,ed,beta(inu),eta(inu),status)
          if (status.ne.0) then
            beta(inu) = 0.
            eta(inu) = 0.
          endif
        enddo

!.....calculate distribution function...................................
        call spectrum_add(i,d(i),ynu(:,i),znu(:,i),iray,localft,statnut)
        do inu=1,nut
          if (statnut(inu).ne.0) localft(:,inu) = 0.
        enddo

!.....convert units.....................................................
        do inu=1,nut
          do ie=1,ne
            localft(ie,inu) = localft(ie,inu)*ps(ie)/d(i)
          enddo
        enddo

!.....calculate rates...................................................
        yedot(len) = -sum(sl(:,len,i))                     !1/s
        yedot(lea) =  sum(sl(:,lea,i))                     !1/s
        edot(len) = -sum(sl(:,len,i)*E)*units%MeV/units%mb !erg/g/s
        edot(lea) = -sum(sl(:,lea,i)*E)*units%MeV/units%mb !erg/g/s

!.....calculate luminosity..............................................
        do inu=1,nut
          tmp = srccp(:,inu,i)*d(i)/units%mb*dV(i)
          LN(inu) = LN(inu) + sum(tmp)             !particles/s
          LE(inu) = LE(inu) + sum(tmp*E)*units%MeV !erg/s
          LNk(:,inu,i) = LNk(:,inu,i) + tmp            !particles/s
          LEk(:,inu,i) = LEk(:,inu,i) + tmp*E*units%MeV!erg/s
        enddo

!.....write profiles of conditions......................................
        tnu = 1./(1.d-98+beta)
        write(1,33) i,errflag(i),
     &    ar(i),d(i),t,ye(i),mu(je),mu(jhat),             !matter
     &    ynu(1:2,i),znu(1:2,i),tnu(1:2),eta(1:2),        !neutrinos
     &    yedot(1:2),edot(1:2),LN(1:2),LE(1:2)            !sources
33      format(2i6,24e12.4)
 
        where (localft.lt.0.1d-98)
          localft = 0.
        end where
      	ft(:,:,i) = localft
      end do

!.....write distribution function.......................................
      do inu=1,nut
        do i=1,n
          write(2,33) i,errflag(i),ft(:,inu,i)
          write(3,33) i,errflag(i),fs(:,inu,i)
          write(4,33) i,errflag(i),lambda(:,inu,i)
        end do	 
      end do

!.....write luminosities................................................
      do inu=1,nut
        do i=1,n
          write(7,33) i,errflag(i),LNk(:,inu,i)
          write(8,33) i,errflag(i),LEk(:,inu,i)
        enddo
      enddo

      close(1)
      close(2)
      close(3)
      close(4)
      close(7)
      close(8)
      deallocate(ar)

      end subroutine wrnuprox

!=======================================================================
