!.....read in ec-rates..................................................
      subroutine read_ec

      use egroup_module, only: ne,e,egroup_empty,egroup
      use ec_module, only: itmax,iymax,idmax,iemax,
     &                     ec_d,ec_t,ec_ye,ec,
     &                     ec_erbmin,ec_erbmax
      use input_file_module
      use units_module

      implicit none

!.....locals............................................................
      integer :: i,k,it,iy,id,ie
      real  :: dea,re,norm,de
      real, dimension(iemax) :: ec_erb
      real, dimension(idmax,iymax,itmax,iemax) :: ec0
      character(len=4) :: duma

!.....build ec-rate tabel v-energy grid.................................
      ec_erb(1) = ec_erbmin
      ec_erb(iemax) = ec_erbmax
      dea = (ec_erb(iemax)-ec_erb(1))/float(iemax-1)
      do i=2,iemax
        ec_erb(i) = ec_erb(i-1) + dea
      enddo
      dea=0.0

!.....open ec-rate table filename.......................................
      if (input_path%initial) call input_read_path
      open(unit=50,file=trim(input_path%ec),status='old')

!.....skip header.......................................................
      do i=1,23
        read(50,*)
      enddo

!.....read table........................................................
      do id=1,idmax
        do iy=1,iymax
          read(50,51) duma,ec_d(id),ec_ye(iy)
          do it=1,itmax
            read(50,*) ec_t(it),norm,ec0(id,iy,it,1:iemax)

!.....set emissivity [1/cm].............................................
            ec0(id,iy,it,1:iemax) = norm*ec0(id,iy,it,1:iemax)
     &                            / (ec_erb(1:iemax)*ec_erb(1:iemax))
     &                            * (units%h*units%c)**3/(2.*units%pi)
     &                            * ec_d(id)/units%mb
     &                            / units%c
         enddo
        enddo
      enddo
51    format(a4,f6.1,5x,f6.2)
!52    format(f6.1,1x,e12.5,1x,201(1x,pe10.4))
      close(50)

!.....build new ec-rate table with IDSA energy bins.....................
      if (egroup_empty) call egroup
      do id=1,idmax
        do iy=1,iymax
          do it=1,itmax
            do k=1,ne

              if (e(k).le.ec_erb(iemax)) then

                if (e(k).le.ec_erb(1)) then
                  ie=1
                  re=0.
                  write(6,85) k,re
                  write(6,86) e(k),ec_erb(1)
                elseif (e(k).ge.ec_erb(iemax)) then
                  ie=iemax-1
                  re=1.
                  write(6,85) k,re
                  write(6,86) e(k),ec_erb(iemax)
                else
                  do i=2,iemax-1
                    if ( (e(k).ge.ec_erb(i)).and.
     &                   (e(k).lt.ec_erb(i+1)) )then
                      ie = i
                      de = ec_erb(ie+1) - ec_erb(ie)
                      re = (e(k) - ec_erb(ie))/de
                    endif
                  enddo
                endif
                ec(id,iy,it,k) = (1.-re)*ec0(id,iy,it,ie) 
     &                         +     re *ec0(id,iy,it,ie+1)
              else
                ec(id,iy,it,k) = 0.
              endif

            enddo            
          enddo            
        enddo            
      enddo            

85    format(' Warning: ec table out of range: re=',f10.4)
86    format(' (ev,ev_ec)',2(2x,e11.4))


!.....debugging.........................................................
      write(6,54) 'T: ',ec_t(1:itmax)
      write(6,54) 'd: ',ec_d(1:idmax)
      write(6,54) 'ye:',ec_ye(1:iymax)
      write(6,54) 'E: ',ec_erb(1:iemax)
      write(6,*)
54    format(a5,210(1x,e12.4))

      write(6,55) 'EC:',ec(1,1,1,1:5)
      write(6,55) 'EC:',ec(1,1,2,1:5)
      write(6,55) 'EC:',ec(1,2,1,1:5)
      write(6,55) 'EC:',ec(idmax,iymax,itmax,1:5)
55    format(a5,30(1x,pe10.4))

      end subroutine



!.....interpolation routine.............................................
      subroutine interpolate_ec(t9,ye,d,em)
 
      use egroup_module, only: ne,e
      use ec_module, only: itmax,iymax,idmax,iemax,
     &                     ec_t,ec_d,ec_ye,ec

      implicit none

!.....input variables data for density, temperature & electron fraction.
      real :: d    ! [g/cm^3]
      real :: t9
      real :: ye

!.....output............................................................
      real, dimension(ne) :: em

!----locals ------------------------------------------------------------
      real, dimension(itmax) :: tt
      real, dimension(iymax) :: ty
      real, dimension(idmax) :: td

      real :: dt,dy,dd,de
      real :: rd,rt,ry,re
      integer :: it,iy,id
      integer :: i,k

!.....find t-index......................................................
      if (t9.le.ec_t(1)) then
        it=1
        rt=0.0
        write(6,81) rt
        write(6,82) t9,ye,d
      elseif (t9.ge.ec_t(itmax)) then
        it=itmax-1
        rt=1.0
        write(6,81) rt
        write(6,82) t9,ye,d
      else
        do i=2,itmax-1
          if ( (t9.ge.ec_t(i)).and.(t9.lt.ec_t(i+1)) )then
            it = i
            dt = ec_t(it+1) - ec_t(it)
            rt = (t9 - ec_t(it))/dt
            exit
          endif
        enddo
      endif
      if (it.lt.1) then
        it=1
        rt=0.
      endif
      if (it.gt.itmax) then
        it=itmax
        rt=0.
      endif

!.....find ye-index.....................................................
      if (ye.le.ec_ye(1)) then
        iy=1
        ry = 0.0
        write(6,83) ry
        write(6,82) t9,ye,d
      elseif (ye.ge.ec_ye(iymax)) then
        iy=iymax-1
        ry=1.0
        write(6,83) ry
        write(6,82) t9,ye,d
      else
        do i=2,iymax-1
          if ( (ye.ge.ec_ye(i)).and.(ye.lt.ec_ye(i+1)) )then
            iy = i
            dy = ec_ye(iy+1) - ec_ye(iy)
            ry = (  ye - ec_ye(iy))/dy
            exit
          endif
        enddo
      endif

!.....find d-index......................................................
      if (d.le.ec_d(1)) then
        id=1
        rd=0.0
        write(6,84) rd
        write(6,82) t9,ye,d
      elseif (d.ge.ec_d(idmax)) then
        id=idmax-1
        rd=1.0
        write(6,84) rd
        write(6,82) t9,ye,d
      else
        do i=1,idmax-1
          if ( (d.ge.ec_d(i)).and.(d.lt.ec_d(i+1)) )then
            id = i
            dd = ec_d(id+1) - ec_d(id)
            rd = (d - ec_d(id))/dd
            exit
          endif
        enddo
      endif

!.TF.debug.
!      write(*,666) t9,ye,d
!      write(*,*) itmax,iymax,idmax
!      write(*,*) it,iy,id
!666   format(i4,3(2x,e12.4))

81    format(' Warning: ec table out of range: rt=',f10.4)
82    format(' (T,ye,rho)',3(2x,e11.4))
83    format(' Warning: ec table out of range: ry=',f10.4)
84    format(' Warning: ec table out of range: rd=',f10.4)


!.....interpolation scheme.............................................
      do k=1,ne
        em(k) =
     &    (1.-rt)*(  (1.-ry)*
     &                       ( (1.-rd)*(ec(id  ,iy  ,it  ,k))
     &   +                         rd *(ec(id+1,iy  ,it  ,k)) )
     &   +               ry *( (1.-rd)*(ec(id  ,iy+1,it  ,k))
     &   +                         rd *(ec(id+1,iy+1,it  ,k)) ) )
     &   +    rt *(  (1.-ry)*( (1.-rd)*(ec(id  ,iy  ,it+1,k))
     &   +                         rd *(ec(id+1,iy  ,it+1,k)) )
     &   +               ry *( (1.-rd)*(ec(id  ,iy+1,it+1,k))
     &   +                         rd *(ec(id+1,iy+1,it+1,k)) ) )
      enddo

      end subroutine interpolate_ec
