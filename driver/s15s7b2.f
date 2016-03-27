!=======================================================================
!
!     Reading progenitor model s15s7b2
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine s15s7b2(am,ar,au,bd,bt,by)

      use agile_parameter_module
      use input_file_module
      use units_module
 
      implicit none

      real, dimension(nnp), intent(out) :: am,ar,au,bd,bt,by

!-----------------------------------------------------------------------

      integer :: i,ip,np
      real :: tmp1,tmp2
      real, dimension(nnp+1) :: bnm,anr,bnt,bnd,anu,bny

!.....initialise input..................................................
      if (input_path%initial) call input_read_path
      if (input_state%initial) call input_read_state
      np = input_state%np

!.....read in from file.................................................
     
      open(unit=10,file=trim(input_path%model)
     &  //trim(input_path%progenitorfn),status='old')


!....s15/s40
!      do ip=1,np+1
!        read(10,11) i,bnm(ip),anr(ip),bnt(ip),bnd(ip),anu(ip),bny(ip)
!      enddo
!11    format(i8,6g12.5)

!....h40
!      open(unit=10,file=trim(input_path%model)
!     &  //trim(input_path%progenitorfn),status='old')
!      do ip=1,np+1
!        read(10,*) i,bnm(ip),anr(ip),bnt(ip),bnd(ip),anu(ip),bny(ip)
!      enddo
!11    format(i8,6g12.5)

!....for u40/u50...............................................................
      read(10,*)
      do ip=1,np+1
        read(10,*) i,bnm(ip),bnd(ip),bnt(ip),anr(ip),bny(ip),anu(ip)
        bnd(ip)=10.0**bnd(ip)
        bnt(ip)=10.0**bnt(ip)
!        write(6,666) i,bnm(ip),bnd(ip),bnt(ip),anr(ip),bny(ip),anu(ip)
      enddo
!666   format(i4,10(2x,e12.5))
      close(10)


!...s15/s40/u40/u50  with units
!...h40 without


!-----shift a-grid down half a zone-------------------------------------
      do ip=1,np
        am(ip) = bnm(ip)*units%Ms                                   ![g]
      enddo

!.....average radius....................................................
      ar(1) = 0.
      tmp1 = anr(1)**3
      do ip=2,np
        tmp2 = anr(ip)**3
        ar(ip) = (0.5*(tmp1 + tmp2))**(1./3.)
        tmp1 = tmp2
      enddo

!.....average velocity..................................................
      au(1) = 0.
      tmp1 = anr(1)**2*anu(1)
      do ip=2,np
        tmp2 = anr(ip)**2*anu(ip)
        au(ip) = 0.5*(tmp1 + tmp2)/ar(ip)**2
        tmp1 = tmp2
      enddo
        
!-----shift b-grid up half a zone---------------------------------------

!.....average density...................................................
      do ip=1,np
        bd(ip) = 10.**(0.5*(log10(bnd(ip+1)) + log10(bnd(ip))))
      enddo
      
!.....average temperature...............................................
      do ip=1,np
        bt(ip) = 10.**(0.5*(log10(bnt(ip+1)) + log10(bnt(ip))))
      enddo

!.....average electron fraction.........................................
      do ip=1,np
        by(ip) = 0.5*(bny(ip+1) + bny(ip))
      enddo

!.....integrate radius from density profile.............................
      ar(1)=0.0
      do ip=2,np
        ar(ip)=(ar(ip-1)**3+(am(ip)-am(ip-1))
     &    /(4.*units%pi*bd(ip-1)/3.))**(1.0/3.0) 
      enddo

      end subroutine s15s7b2

!=======================================================================
