!=======================================================================
!
!     AGILE: atvdflux
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine btvdflux(nq,ada,bda,bv,aw,flx)

      use agile_parameter_module

      implicit none

      integer, intent(in) :: nq
      real, dimension(nnq), intent(in) :: ada,bda,bv,aw
      real, dimension(nnq), intent(out) :: flx

!-----------------------------------------------------------------------

      integer :: iq
      real :: w1,w2,w3,df1p,df2p,d1,d2,tmp,v1p,da1p,da2p

!.....initialisation....................................................
      w2 = 0.
      w3 = 0.
      df2p = 0.
      d2 = 0.
      da2p = 0.
      do iq=1,nq

!.....flux..............................................................
        w1 = w2
        w2 = w3
        w3 = aw(iq)
        if (iq.eq.1) cycle

!.....flux correction...................................................
        da1p = da2p
        da2p = bda(iq)
        df1p = df2p
        df2p = (w3-w2)/da2p
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
        v1p = bv(iq-1)
        if (v1p.gt.0.) then
          flx(iq-1) = (w1 + d1*0.5*da1p)*v1p
        else
          flx(iq-1) = (w2 - d2*0.5*da1p)*v1p
        endif
      enddo

!.....boundary..........................................................
      flx(1) = bv(1)*aw(1)
      flx(nq) = bv(nq)*aw(nq)

      end subroutine btvdflux

!=======================================================================
