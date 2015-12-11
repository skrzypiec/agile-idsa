!=======================================================================
!
!     AGILE: Equations for relax on initial profile
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine equations_relax(jacobi,o,n,f,status)

      use adaptive_grid_module
      use agile_parameter_module
      use state_vector_module

      logical, intent(in) :: jacobi
      integer, intent(out) :: status
      real, dimension(100) :: ut
      real, dimension(nnc), intent(out) :: f
      type(state_vector), intent(in) :: o,n

!-----------------------------------------------------------------------
!
!     input:
!     jacobi   is true if the routine is called for derivatives
!     o   old state vector (used for boundaries)
!     n   guess state vector
!
!     output:
!     f   residua
!     status is set to 1 in case of a problem
!
!     This routine is written for voffset=4.
!     There are ny-voffset boundary conditions at left and voffset
!     boundary conditions at right hand. The grid equations is the
!     (voffset+1)-th equation in the view of a zone.
!
!-----------------------------------------------------------------------

      logical, save :: initial=.true.
      integer, save :: nq,nq1,nz,nz1
      integer :: ie,i
      integer, dimension(nnq+1,nny), save :: jf
      real, dimension(nnq) :: atmp,btmp
      type(state_vector) :: ref

!-----initialise--------------------------------------------------------

!.....abbreviations for array-lengths...................................
      if (initial) then
        nq = n%nq
        nq1 = nq-1
        nz = nq+1
        nz1 = nz-1

!.....assign zones to residuum vector...................................
        if (input_jacobi%initial) call input_read_jacobi
        do ie=input_jacobi%voffset+1,n%ny
          jf(1:nz1,ie) = (/(i,i=0,nz-2)/)*n%ny+ie-input_jacobi%voffset
        enddo
        do ie=1,input_jacobi%voffset
          jf(2:nz,ie) = (/(i,i=1,nz1)/)*n%ny+ie-input_jacobi%voffset
        enddo
        initial = .false.
      endif
      
!-----precalculate some terms for later use-----------------------------
    
!.....create reference state ref from profile...........................
      ref = n
      call adaptive_linint(ref,status)
      if (status.ne.0) return
      call state_fill_zones(ref,status)
      if (status.ne.0) return
      
!-----calculate equations-----------------------------------------------

!.....center boundary...................................................
      f(jf(1,9)) = n%ax(1)
      f(jf(1,10)) = n%ar(1) - o%ar(1)
      f(jf(1,11)) = n%au(1) - o%au(1)
      f(jf(1,12)) = n%am(1) - o%am(1)

!.....center ghost zones................................................
      f(jf(2,1)) = n%bin%n(2) - n%bin%n(1)
      f(jf(2,2)) = n%bin%t(2) - n%bin%t(1)
      f(jf(2,4)) = n%bin%ye(2) - n%bin%ye(1)
      f(jf(2,3)) = n%balpha(2) - n%balpha(1)
      f(jf(2,5:6)) = n%bynu(2,1:2) - n%bynu(1,1:2)
      f(jf(2,7:8)) = n%bznu(2,1:2) - n%bznu(1,1:2)

!.....interpolate density from profile..................................
      f(jf(3:nz1,1)) = n%bin%n(2:nq1) - ref%bin%n(2:nq1)

!.....interpolate temperature from profile..............................
      f(jf(3:nz1,2)) = n%bin%t(2:nq1) - ref%bin%t(2:nq1)

!.....interpolate velocity from profile.................................
      f(jf(3:nz1,3)) = n%au(2:nq1) - ref%au(2:nq1)
      
!.....interpolate electron fraction from profile........................
      f(jf(3:nz1,4)) = n%bin%ye(2:nq1) - ref%bin%ye(2:nq1)

!.....interpolate neutrino fractions from profile.......................
      f(jf(3:nz1,5)) = n%bynu(2:nq1,1) - ref%bynu(2:nq1,1)
      f(jf(3:nz1,6)) = n%bynu(2:nq1,2) - ref%bynu(2:nq1,2)
      f(jf(3:nz1,7)) = n%bznu(2:nq1,1) - ref%bznu(2:nq1,1)
      f(jf(3:nz1,8)) = n%bznu(2:nq1,2) - ref%bznu(2:nq1,2)

!.....volume integration................................................
      btmp(2:nq) = n%aV(2:nq) - n%aV(1:nq1)
      f(jf(2:nz1,10)) = n%bda(2:nq) - n%bD(2:nq)*btmp(2:nq)

!.....gravitational mass integration....................................
      btmp(2:nq) = n%am(2:nq) - n%am(1:nq1)
      if (input_action%rel) then
        f(jf(2:nz1,11))=btmp(2:nq)-(n%bg(2:nq)+n%bei(2:nq))*n%bda(2:nq)
      else
        f(jf(2:nz1,11)) = btmp(2:nq)-n%bda(2:nq)
      endif

!.....lapse function....................................................
      if (input_action%rel) then
        btmp(1:nq) = n%bD(1:nq)*(n%bg(1:nq) + n%bei(1:nq))     ![g/cm^3]
        atmp(2:nq1) = 0.5*(btmp(3:nq) + btmp(2:nq1))           ![g/cm^3]
        f(jf(2:nz-2,12)) = (n%bp(3:nq) - n%bp(2:nq1))/units%c**2
     1    + atmp(2:nq1) * (n%balpha(3:nq) - n%balpha(2:nq1))   ![g/cm^3]
      else
        f(jf(2:nz-2,12)) = n%balpha(3:nq) - n%balpha(2:nq1)          ![]
      endif

!.....surface ghost.....................................................
      f(jf(nz,1)) = n%bin%n(nq) - o%bin%n(nq)
      f(jf(nz,2)) = n%bin%t(nq) - o%bin%t(nq)
      f(jf(nz,3)) = n%balpha(nq)*n%aw(nq)
     1            - 0.5*(n%balpha(nq)+n%balpha(nq1))*n%aw(nq1)
      f(jf(nz,4)) = n%bin%ye(nq) - o%bin%ye(nq)
      f(jf(nz,5:6)) = n%bynu(nq,1:2) - o%bynu(nq,1:2)
      f(jf(nz,7:8)) = n%bznu(nq,1:2) - o%bznu(nq,1:2)

!.....surface boundary..................................................
      f(jf(nz1,9)) = n%ax(nq)
      if (input_action%rel) then
        f(jf(nz1,12)) = n%balpha(nq) - (1.-2.*n%amr(nq))/n%ag(nq)
      else
        f(jf(nz1,12)) = n%balpha(nq) - 1.
      endif

      end subroutine equations_relax

!=======================================================================
