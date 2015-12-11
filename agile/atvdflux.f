!=======================================================================
!
!     AGILE: atvdflux
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine atvdflux(nq,ada,bda,av,bw,flx)

      use agile_parameter_module

      implicit none

      integer, intent(in) :: nq
      real, dimension(nnq), intent(in) :: ada,bda,av,bw
      real, dimension(nnq), intent(out) :: flx

!-----------------------------------------------------------------------

      integer :: iq
      real :: w1,w2,w3,df1p,df2p,d1,d2,tmp,v1p,da1,da2

!.....initialisation....................................................
      w2 = 0.
      w3 = 0.
      df2p = 0.
      d2 = 0.
      da2 = 0.
      do iq=1,nq

!.....flux..............................................................
        w1 = w2
        w2 = w3
        w3 = bw(iq)
        if (iq.eq.1) cycle

!.....flux correction...................................................
        df1p = df2p
        df2p = (w3-w2)/ada(iq-1)
        if (iq.eq.2) cycle

!.....Van Leer..........................................................
        d1 = d2
        tmp = df2p*df1p
        if (tmp.gt.0.) then
          d2 = 2.*tmp/(df2p+df1p)
        else
          d2 = 0.
        endif

!.....upwind............................................................
        v1p = av(iq-2)
        da1 = da2
        da2 = bda(iq-1)
        if (v1p.gt.0.) then
          flx(iq-2) = (w1 + d1*0.5*da1)*v1p
        else
          flx(iq-2) = (w2 - d2*0.5*da2)*v1p
        endif
      enddo

!.....zone nq-1.........................................................
      w1 = w2
      w2 = w3
      v1p = av(nq-1)
      if (v1p.gt.0.) then
        flx(nq-1) = w1*v1p
      else
        flx(nq-1) = w2*v1p
      endif

!.....boundary..........................................................
      flx(1) = av(1)*bw(1)
      flx(nq) = av(nq)*bw(nq)

      end subroutine atvdflux

!=======================================================================
