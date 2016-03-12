!***********************************************************************
!
!     eos_module maintains a lookup table for the equation of state
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module eos_module

      use units_module

      implicit none

!=======================================================================
!
!     EOS: parameters
!
!=======================================================================

      integer, parameter :: nnx=3       !number of input quantities
      integer, parameter :: nn=130      !number of table entries
      integer, parameter :: nny=11      !number of output quantities

!=======================================================================
!
!     EOS: type definitions
!
!=======================================================================

      type eos_table
        logical :: initial               !flag if table is loaded
        logical, dimension(nnx) :: lgx   !choice of interpolation method
        logical, dimension(nny) :: lgy   !choice of interpolation method
        real(4), dimension(nnx) :: xmin,dx  !table min and increment
        real(4), dimension(nn,nn,nn) :: e   !internal energy table
        integer, dimension(nn,nn,nn) :: ef !eos flag
        real(4), dimension(nn,nnx) :: xtab  !grid vector
        real(4), dimension(nny,nn,nn,nn) :: y !table with eos output
      end type eos_table

      type eos_state
        real :: t,e,s,p,a,z,xa,xp,xn,muhat,mun,mue
      end type eos_state

!=======================================================================
!
!     EOS: eos table
!
!=======================================================================

      real, parameter :: ezero=8.9 !zero energy in MeV/(baryon mass)
      real, parameter :: MASSN=939.5731 !neutron mass used in LS-EoS
      real, parameter :: eionglobal=(MASSN-ezero)*units%MeV/units%mb
     &  - units%c**2                                         ![cm^2/s^2]
      type(eos_table) :: ls

      contains

!=======================================================================
!
!     EOS: read table from file
!
!=======================================================================

      subroutine eos_read(filename)

      implicit none

      character(*), intent(in) :: filename

!-----------------------------------------------------------------------
!
!     Note: The definition of the internal energy depends on the
!     definition of the baryon mass. If the eos table is the
!     Lattimer-Swesty equation of state, the 'internal' energy in the
!     table is (where UTOT and MASSN as defined as in the LS eos)
!     e = (UTOT+MASSN)*units%MeV/units%mb - units%c**2 [erg/g].
!
!-----------------------------------------------------------------------

      logical :: found
      integer :: is,iy,id,istart,count

      inquire(file=trim(filename),exist=found)
      if (.not.found) write(6,*) 'Error: File ',trim(filename),
     &  ' not found!'
      open(1,file=trim(filename),form='unformatted',status='old')
      read(1) ls
      close(1)
      ls%initial = .false.

!.....eliminating holes altogether......................................
      count = 0
      do id=1,nn
        do iy=1,nn
          is = 1
          do while (ls%ef(is,iy,id).eq.0)
            is = is+1
            if (is.gt.nn) then
              write(6,*) 'Warning: cannot fill hole in eos_module.f'
              write(6,*) 'id,iy: ',id,iy
              exit
            endif
          enddo
          istart = is
          if (is.gt.1) count = count+1
          do is=istart+1,nn
            if (ls%ef(is,iy,id).eq.0) then
              ls%e(is,iy,id) = ls%e(is-1,iy,id)
              ls%y(:,is,iy,id) = ls%y(:,is-1,iy,id)
              ls%ef(is,iy,id) = ls%ef(is-1,iy,id)
              if (ls%ef(is,iy,id).eq.0)
     &          write(6,*) 'Warning: Unfixed EoS hole at: ',is,iy,id
            endif
          enddo
          do is=istart-1,1,-1
            if (ls%ef(is,iy,id).eq.0) then
              ls%e(is,iy,id) = ls%e(is+1,iy,id)
              ls%y(:,is,iy,id) = ls%y(:,is+1,iy,id)
              ls%ef(is,iy,id) = ls%ef(is+1,iy,id)
              if (ls%ef(is,iy,id).eq.0)
     &          write(6,*) 'Warning: EoS hole at: ',is,iy,id
            endif
          enddo
        enddo
      enddo
      if (count.gt.0) write(6,11) count
11    format('Warning: ',i9,' holes fixed in EoS table')

      write(6,60) ' Info: d-range  =',ls%xmin(1),ls%dx(1)
      write(6,60) ' Info: ye-range =',ls%xmin(2),ls%dx(2)
      write(6,60) ' Info: s-range  =',ls%xmin(3),ls%dx(3)
60    format(a16,2x,3(f13.5))


      end subroutine eos_read

!=======================================================================
!
!     EOS: interpolate in eos table with energy input
!
!=======================================================================

      subroutine eos_einterp(d,ye,e,cs,state,status)

      implicit none

      integer, intent(out) :: status
      real, intent(in) :: d,ye,e
      real, intent(out) :: cs
      type(eos_state), intent(out) :: state

!-----------------------------------------------------------------------
!
!     input:
!     d ... density [g/cm^3]
!     ye ... electron fraction []
!     e ... specific internal energy [erg/g]
!
!     output:
!     s ... entropy [kB/baryon]
!     p ... pressure [erg/cm^3]
!     cs ... sound speed [cm/s]
!     t ... temperature [MeV]
!     status ... >0 on failure
!
!     Note: The definition of the internal energy depends on the 
!     definition of the baryon mass. If the eos table is the
!     Lattimer-Swesty equation of state, the 'internal' energy in the
!     table is (where UTOT and MASSN as defined as in the LS eos)
!     e = (UTOT+MASSN)*units%MeV/units%mb - units%c**2 [erg/g].
!
!-----------------------------------------------------------------------

      integer :: n,il,iu,im,id,iy,is
      integer, dimension(3) :: ix
      real :: logd,rd,ry,rs
      real, dimension(nn) :: etmp
      real, dimension(nny) :: y

      status = 0
      n = nn
      if (ls%initial) then
        write(6,*) 'Error: eos table unknown in eos_einterp!'
        stop
      endif

!.....find density index................................................
      logd = log10(d)
      id = 1 + floor( (logd - ls%xmin(1))/ls%dx(1) )
      if (id.lt.1) then
        status = 1
        id = 1
        rd = 0.
      elseif (id.ge.nn) then
        status = 1
        id = nn-1
        rd = 1.
      else
        rd = (logd - ls%xtab(id,1))/ls%dx(1)
      endif

!.....find ye index.....................................................
      iy = 1 + floor( (ye - ls%xmin(2))/ls%dx(2) )
      if (iy.lt.1) then
        status = 2
        iy = 1
        ry = 0.
      elseif (iy.ge.nn) then
        status = 2
        iy = nn-1
        ry = 1.
      else
        ry = (ye - ls%xtab(iy,2))/ls%dx(2)
      endif

!.....find entropy index................................................
      il = 0
      iu = n+1
      do while (iu-il.gt.1)
        is = (iu+il)/2
        etmp(is) = (1.-rd)*( (1.-ry)*ls%e(is  ,iy  ,id  )
     1           +               ry *ls%e(is  ,iy+1,id  ) )
     2           +     rd *( (1.-ry)*ls%e(is  ,iy  ,id+1)
     3           +               ry *ls%e(is  ,iy+1,id+1) )
        if (e.gt.etmp(is)) then
          il = is
        else
          iu = is
        endif
      enddo
      is = il
      if (is.lt.1) then
        status = 3
        is = 1
        rs = 0.
      elseif (is.ge.nn) then
        status = 3
        is = nn-1
        rs = 1.
      else
        rs = (e - etmp(is))/(etmp(is+1) - etmp(is))
      endif

!.....check for holes in eos table......................................
      if (any(ls%ef(is:is+1,iy:iy+1,id:id+1).eq.0)) then
        write(6,*) 'Warning: There is a hole in the EoS table!'
        write(6,*) 'is,iy,id: ',is,iy,id
        status = 4
        return
      endif

!.....interpolate.......................................................
      y = (1.-rd)*( (1.-ry)*( (1.-rs)*ls%y(1:nny,is  ,iy  ,id  )
     1  +                         rs *ls%y(1:nny,is+1,iy  ,id  ) )
     2  +               ry *( (1.-rs)*ls%y(1:nny,is  ,iy+1,id  )
     3  +                         rs *ls%y(1:nny,is+1,iy+1,id  ) ) )
     4  +     rd *( (1.-ry)*( (1.-rs)*ls%y(1:nny,is  ,iy  ,id+1)
     5  +                         rs *ls%y(1:nny,is+1,iy  ,id+1) )
     6  +               ry *( (1.-rs)*ls%y(1:nny,is  ,iy+1,id+1)
     7  +                         rs *ls%y(1:nny,is+1,iy+1,id+1) ) )

!.....assign results and return.........................................
      state%e = e
      state%p = 10.**y(3)
      cs = sqrt(abs(ls%y(3,is+1,iy+1,id+1)-ls%y(3,is+1,iy+1,id))
     &   / ls%dx(1)*state%p/d)
      state%t = y(2)
      state%s = 10.**(ls%xtab(is,3)+rs*ls%dx(3))
      state%a = y(4)
      state%z = y(5)
      state%xa = y(6)
      state%xp = y(7)
      state%xn = y(8)
      state%muhat = y(9)
      state%mun = y(10)
      state%mue = y(11)

      end subroutine eos_einterp

!=======================================================================
!
!     EOS: interpolate in eos table with temperature input
!
!=======================================================================

      subroutine eos_tinterp(d,ye,t,cs,state,status)

      implicit none

      integer, intent(out) :: status
      real, intent(in) :: d,ye,t
      real, intent(out) :: cs
      type(eos_state), intent(out) :: state

!-----------------------------------------------------------------------
!
!     input:
!     d ... density [g/cm^3]
!     ye ... electron fraction []
!     t ... temperature [MeV]
!
!     output:
!     s ... entropy [kB/baryon]
!     p ... pressure [erg/cm^3]
!     cs ... sound speed [cm/s]
!     e ... specific internal energy [erg/g]
!     status ... >0 on failure
!
!-----------------------------------------------------------------------

      integer :: n,il,iu,im,id,iy,is
      integer, dimension(3) :: ix
      real :: logd,rd,ry,rs
      real, dimension(nn) :: ttmp
      real, dimension(nny) :: y

      status = 0
      n = nn
      if (ls%initial) then
        write(6,*) 'Error: eos table unknown in eos_tinterp!'
        stop
      endif

!.....find density index................................................
      logd = log10(d)
      id = 1 + floor( (logd - ls%xmin(1))/ls%dx(1) )
      if (id.lt.1) then
        status = 1
        id = 1
        rd = 0.
      elseif (id.ge.nn) then
        status = 1
        id = nn-1
        rd = 1.
      else
        rd = (logd - ls%xtab(id,1))/ls%dx(1)
      endif

!.....find ye index.....................................................
      iy = 1 + floor( (ye - ls%xmin(2))/ls%dx(2) )
      if (iy.lt.1) then
        status = 2
        iy = 1
        ry = 0.
      elseif (iy.ge.nn) then
        status = 2
        iy = nn-1
        ry = 1.
      else
        ry = (ye - ls%xtab(iy,2))/ls%dx(2)
      endif

!.....find temperature index............................................
      il = 0
      iu = n+1
      do while (iu-il.gt.1)
        is = (iu+il)/2
        ttmp(is) = (1.-rd)*( (1.-ry)*ls%y(2,is  ,iy  ,id  )
     1           +               ry *ls%y(2,is  ,iy+1,id  ) )
     2           +     rd *( (1.-ry)*ls%y(2,is  ,iy  ,id+1)
     3           +               ry *ls%y(2,is  ,iy+1,id+1) )
        if (t.gt.ttmp(is)) then
          il = is
        else
          iu = is
        endif
      enddo
      is = il
      if (is.lt.1) then
        status = 3
        is = 1
        rs = 0.
      elseif (is.ge.nn) then
!.....TF: neglect eos-table T-index exceed for b.h. formation...........
!        status = 3
        is = nn-1
        rs = 1.
      else
        rs = (t - ttmp(is))/(ttmp(is+1) - ttmp(is))
      endif

!.....check for holes in eos table......................................
      if (any(ls%ef(is:is+1,iy:iy+1,id:id+1).eq.0)) then
        write(6,*) 'Warning: There is a hole in the EoS table!'
        write(6,*) 'is,iy,id: ',is,iy,id
        status = 4
        return
      endif

!.....interpolate.......................................................
      y = (1.-rd)*( (1.-ry)*( (1.-rs)*ls%y(1:nny,is  ,iy  ,id  )
     1  +                         rs *ls%y(1:nny,is+1,iy  ,id  ) )
     2  +               ry *( (1.-rs)*ls%y(1:nny,is  ,iy+1,id  )
     3  +                         rs *ls%y(1:nny,is+1,iy+1,id  ) ) )
     4  +     rd *( (1.-ry)*( (1.-rs)*ls%y(1:nny,is  ,iy  ,id+1)
     5  +                         rs *ls%y(1:nny,is+1,iy  ,id+1) )
     6  +               ry *( (1.-rs)*ls%y(1:nny,is  ,iy+1,id+1)
     7  +                         rs *ls%y(1:nny,is+1,iy+1,id+1) ) )

!.....assign results and return.........................................
      state%t = t
      state%p = 10.**y(3)
      cs = sqrt(abs(ls%y(3,is+1,iy+1,id+1)-ls%y(3,is+1,iy+1,id))
     &   / ls%dx(1)*state%p/d)
      state%e = y(1)
      state%s = 10.**(ls%xtab(is,3)+rs*(ls%xtab(is+1,3)-ls%xtab(is,3)))
      state%a = y(4)
      state%z = y(5)
      state%xa = y(6)
      state%xp = y(7)
      state%xn = y(8)
      state%muhat = y(9)
      state%mun = y(10)
      state%mue = y(11)

      end subroutine eos_tinterp

!=======================================================================
!
!     EOS: interpolate in eos table with entropy input
!
!=======================================================================

      subroutine eos_sinterp(d,ye,s,cs,state,status)

      implicit none

      integer, intent(out) :: status
      real, intent(in) :: d,ye,s
      real, intent(out) :: cs
      type(eos_state), intent(out) :: state

!-----------------------------------------------------------------------
!
!     input:
!     d ... density [g/cm^3]
!     ye ... electron fraction []
!     s ... entropy [kB/baryon]
!
!     output:
!     e ... specific internal energy [erg/g]
!     p ... pressure [erg/cm^3]
!     cs ... sound speed [cm/s]
!     t ... temperature [MeV]
!     status ... >0 on failure
!
!     Note: The definition of the internal energy depends on the 
!     definition of the baryon mass. If the eos table is the
!     Lattimer-Swesty equation of state, the 'internal' energy in the
!     table is (where UTOT and MASSN as defined as in the LS eos)
!     e = (UTOT+MASSN)*units%MeV/units%mb - units%c**2 [erg/g].
!
!-----------------------------------------------------------------------

      integer :: n,il,iu,im,id,iy,is
      integer, dimension(3) :: ix
      real :: logd,logs,rd,ry,rs
      real, dimension(nn) :: etmp
      real, dimension(nny) :: y

      status = 0
      n = nn
      if (ls%initial) then
        write(6,*) 'Error: eos table unknown in eos_sinterp!'
        stop
      endif

!.....find density index................................................
      logd = log10(d)
      id = 1 + floor( (logd - ls%xmin(1))/ls%dx(1) )
      if (id.lt.1) then
        status = 1
        id = 1
        rd = 0.
      elseif (id.ge.nn) then
        status = 1
        id = nn-1
        rd = 1.
      else
        rd = (logd - ls%xtab(id,1))/ls%dx(1)
      endif

!.....find ye index.....................................................
      iy = 1 + floor( (ye - ls%xmin(2))/ls%dx(2) )
      if (iy.lt.1) then
!        status = 2
        iy = 1
        ry = 0.
      elseif (iy.ge.nn) then
!        status = 2
        iy = nn-1
        ry = 1.
      else
        ry = (ye - ls%xtab(iy,2))/ls%dx(2)
      endif

!.....find density and ye index.........................................
      logs = log10(s)
      is = 1 + floor( (logs - ls%xmin(3))/ls%dx(3) )
      if (is.lt.1) then
!.....TF: neglect eos limit for black-hole formation....................
!        status = 3
        is = 1
        rs = 0.
      elseif (is.ge.nn) then
!.....TF: neglect eos limit for black-hole formation....................
!        status = 3
        is = nn-1
        rs = 1.
      else
        rs = (logs - ls%xtab(is,3))/ls%dx(3)
      endif

!.....check for holes in eos table......................................
      if (any(ls%ef(is:is+1,iy:iy+1,id:id+1).eq.0)) then
        write(6,*) 'Warning: There is a hole in the EoS table!'
        write(6,*) 'is,iy,id: ',is,iy,id
        status = 4
        return
      endif

!.....interpolate.......................................................
      y = (1.-rd)*( (1.-ry)*( (1.-rs)*ls%y(1:nny,is  ,iy  ,id  )
     1  +                         rs *ls%y(1:nny,is+1,iy  ,id  ) )
     2  +               ry *( (1.-rs)*ls%y(1:nny,is  ,iy+1,id  )
     3  +                         rs *ls%y(1:nny,is+1,iy+1,id  ) ) )
     4  +     rd *( (1.-ry)*( (1.-rs)*ls%y(1:nny,is  ,iy  ,id+1)
     5  +                         rs *ls%y(1:nny,is+1,iy  ,id+1) )
     6  +               ry *( (1.-rs)*ls%y(1:nny,is  ,iy+1,id+1)
     7  +                         rs *ls%y(1:nny,is+1,iy+1,id+1) ) )

!.....assign results and return.........................................
      state%s = s
      state%p = 10.**y(3)
      cs = sqrt(abs(ls%y(3,is+1,iy+1,id+1)-ls%y(3,is+1,iy+1,id))
     &   / ls%dx(1)*state%p/d)
      state%t = y(2)
      state%e = y(1)
      state%a = y(4)
      state%z = y(5)
      state%xa = y(6)
      state%xp = y(7)
      state%xn = y(8)
      state%muhat = y(9)
      state%mun = y(10)
      state%mue = y(11)

      end subroutine eos_sinterp

!=======================================================================

      end module eos_module

!***********************************************************************
