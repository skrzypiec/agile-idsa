!***********************************************************************
!
!     IDSA: frameint_module, interpolates in frame
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module frameint_module

      interface frameint
        module procedure uframeint
        module procedure nframeint
      end interface

      contains

!=======================================================================
!
!     IDSA: uframeint
!
!=======================================================================

      subroutine uframeint(slog,y,sg,yg,xg,x,dxds,dxdy)

      implicit none

      real, intent(in) :: slog,y
      real, dimension(2), intent(in) :: sg,yg
      real, dimension(2,2), intent(in) :: xg
      real, intent(out) :: x,dxds,dxdy

!-----------------------------------------------------------------------

      real :: c1,c2,tmp

      c1 = (slog-sg(1))/(sg(2)-sg(1))
      c2 = (y-yg(1))/(yg(2)-yg(1))

      tmp = xg(2,2)-xg(1,2)-xg(2,1)+xg(1,1)
      dxds = (xg(2,1)-xg(1,1)+c2*tmp)/(sg(2)-sg(1))
      dxdy = (xg(1,2)-xg(1,1)+c1*tmp)/(yg(2)-yg(1))
      x = xg(1,1)+c1*(xg(2,1)-xg(1,1))+c2*(xg(1,2)-xg(1,1))+c1*c2*tmp

      end subroutine uframeint

!=======================================================================
!
!     IDSA: nframeint
!
!=======================================================================

      subroutine nframeint(slog,y,sg,yg,xg,x,dxds,dxdy)

      use egroup_module
      use species_module

      implicit none

      real, intent(in) :: slog,y
      real, dimension(2), intent(in) :: sg,yg
      real, dimension(ne,nut,2,2), intent(in) :: xg
      real, dimension(ne,nut), intent(out) :: x,dxds,dxdy

!-----------------------------------------------------------------------

      real :: c1,c2
      real, dimension(ne,nut) :: tmp

      c1 = (slog-sg(1))/(sg(2)-sg(1))
      c2 = (y-yg(1))/(yg(2)-yg(1))

      tmp = xg(:,:,2,2)-xg(:,:,1,2)-xg(:,:,2,1)+xg(:,:,1,1)
      dxds = (xg(:,:,2,1)-xg(:,:,1,1)+c2*tmp)/(sg(2)-sg(1))
      dxdy = (xg(:,:,1,2)-xg(:,:,1,1)+c1*tmp)/(yg(2)-yg(1))
      x = xg(:,:,1,1)+c1*(xg(:,:,2,1)-xg(:,:,1,1))
     &  + c2*(xg(:,:,1,2)-xg(:,:,1,1))+c1*c2*tmp

      end subroutine nframeint

!=======================================================================

      end module frameint_module

!***********************************************************************
