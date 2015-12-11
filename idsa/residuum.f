!=======================================================================
!
!     IDSA: residuum, calculates residuum and derivatives for Newton-Raphson
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine residuum(dt,u0,y0,sg,yg,ug,slg,x,res,jac)

      use egroup_module
      use frameint_module
      use species_module

      implicit none

      real, intent(in) :: dt,u0,y0
      real, dimension(2), intent(in) :: sg,yg,x
      real, dimension(2,2), intent(in) :: ug
      real, dimension(ne,nut,2,2), intent(in) :: slg
      real, dimension(2), intent(out) :: res
      real, dimension(2,2), intent(out) :: jac

!-----------------------------------------------------------------------
!
!     ML, CITA (2004)
!     dt ... time step [s]
!     u0 ... old specific energy [erg/g]
!     y0 ... old electron fraction [particles/baryon]
!     sg ... log entropy on frame corners [kB/baryon]
!     yg ... electron fraction on frame corners [electrons/baryon]
!     ug ... specific internal energy on frame corners [MeV/baryon]
!     slg ... neutrino source [particles/baryon/s]
!     x ... guess for unknowns (log temperature,electron fraction)
!     res ... residuum (energy update,electron fraction update)
!     jac ... derivative of residuum with respect to x
!
!-----------------------------------------------------------------------

      real :: u,uds,udy,y
      real, dimension(ne,nut) :: sl,dnds,dndy

!.....interpolate.......................................................
      call frameint(x(1),x(2),sg,yg,ug,u,uds,udy)
      call frameint(x(1),x(2),sg,yg,slg,sl,dnds,dndy)
      y = x(2)
      
!.....calculate residuum................................................
      res(1) = (u-u0)/dt + sum((sl(:,len)+sl(:,lea))*E)
      res(2) = (y-y0)/dt + sum(sl(:,len)-sl(:,lea))
      
!.....calculate derivatives.............................................
      jac(1,1) = uds/dt + sum((dnds(:,len)+dnds(:,lea))*E)
      jac(1,2) = udy/dt + sum((dndy(:,len)+dndy(:,lea))*E)
      jac(2,1) =        + sum(dnds(:,len)-dnds(:,lea))
      jac(2,2) = 1./dt  + sum(dndy(:,len)-dndy(:,lea))
       
      end subroutine residuum

!=======================================================================
