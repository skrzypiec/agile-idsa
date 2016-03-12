!=====================================================================
!
!     signature: outputs neutrino signatures
!
!=====================================================================

      program signature

      implicit none

      integer :: c1,c2,cycle,imax,i,dumi
      integer, parameter :: ip=104
      real, dimension(ip) :: r,Nve,Nveb,Lve,Lveb,Eve,Eveb
      real :: lambda,rad,bounce,t,E,EB,NE,NEB,LE,LEB,dum
      character(len=17) duma
      character(len=25) :: filename
      character(len=13) :: input_filename

!.....filename........................................................
      write(6,*) 'filename?'
      read(5,*) input_filename
      input_filename = 'data/'//input_filename
      write(6,*) 'start/end cycle'
      read(5,*) c1,c2
      write(6,*) 'sampling radius [km]'
      read(5,*) rad
      write(6,*) 'time at bounce [s]'
      read(5,*) bounce

!.....initialization..................................................
      imax = 102
 
!.....open file for output............................................
      open(1,file='neutrino.dat',status='unknown')
      write(1,88)
      write(1,44) rad
44    format('radius R = ',f10.5,' km')
77    format('t',t17,'Le',t47,'Lebar',t47)
      write(1,88)
88    format(116('-'))

      do cycle=c1,c2

!.....read state vector...............................................
        write(filename,'(a13,i5.5,a2)') 'data/nuprox1d',cycle,'.d'
        open(2,file=filename,status='unknown',action='read')
        read(2,*)
        read(2,*)
        read(2,*)
        do i=1,imax
          read(2,*) dumi,dumi,r(i),
     &              dum,dum,dum,dum,dum,dum,dum,dum,dum, 
     &              dum,dum,dum,dum,dum,dum,dum,dum,
     &              Nve(i),Nveb(i),Lve(i),Lveb(i)
          Eve(i)=Lve(i)/Nve(i)*6.2414e5
          Eveb(i)=Lveb(i)/Nveb(i)*6.2414e5
        enddo
        close(2)

        write(filename,'(a9,i4.4,a4)') input_filename,cycle,'.stp'
        open(2,file=filename,status='unknown',action='read')
        read(2,*)
        read(2,*)
        read(2,*)
        read(2,*) duma,t
        close(2)

!.....select radius.....................................................
        i = imax-1
        do while (r(i).gt.rad*1.0e5)
          i = i-1
        enddo
        lambda = ((rad*1.0d5)**3 - r(i)**3)/(r(i+1)**3 - r(i)**3)
        LE  = Lve(i)  + lambda*(Lve(i+1) - Lve(i))
        LEB = Lveb(i) + lambda*(Lveb(i+1) - Lveb(i))
        NE  = Nve(i)  + lambda*(Nve(i+1) - Nve(i))
        NEB = Nveb(i) + lambda*(Nveb(i+1) - Nveb(i))
        E   = Eve(i)  + lambda*(Eve(i+1) - Eve(i))
        EB  = Eveb(i) + lambda*(Eveb(i+1) - Eveb(i))
        write(1,33) t-bounce,LE,LEB,NE,NEB,E,EB
      enddo

33    format(e22.15,6(2x,e15.7))
      close(1)

      end  

!=======================================================================
