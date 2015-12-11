!***********************************************************************
!
!     AGILE: adaptive_grid_module
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module adaptive_grid_module

      use agile_parameter_module
      use input_file_module
      use state_vector_module

      implicit none

!.....here we save the initial profiles for interpolations in linint....
      type(state_vector) :: adaptive_linint_template

      contains

!=======================================================================
!
!     AGILE: equations for relax on allocation
!
!=======================================================================

        subroutine adaptive_relax(jacobi,old,new,f,status)

        logical, intent(in) :: jacobi
        integer, intent(out) :: status
        real, dimension(nnc), intent(out) :: f
        type(state_vector), intent(in) :: old,new

!-----------------------------------------------------------------------
!
!     input:
!     jacobi   is true if the routine is called for derivatives
!     old   old state vector (used for boundaries)
!     new   guess state vector
!
!     output:
!     f   residua
!     status is set to 1 in case of a problem (not used)
!
!-----------------------------------------------------------------------

        integer :: nq,ny,nc,i,ic
        real, dimension(nnq) :: xold,xnew
        type(state_vector) :: tmp

!.....initialise........................................................
        status = 0
        nq = new%nq
        ny = new%ny
        nc = nq*ny

!.....interpolate variables to new grid point locations.................
        tmp = new
        call adaptive_linint(tmp,status)
        if (status.ne.0) return

!.....boundaries of independent variable................................
        call state_sumdx(old,xold)
        call state_sumdx(new,xnew)
        f(1) = xnew(1) - xold(1)
        f(nc-ny+1) = xnew(nq) - xold(nq)
        
!.....other variables...................................................
        do ic=1,nc
          if (mod(ic,ny).ne.1) f(ic) = new%y(ic) - tmp%y(ic)
        enddo

        end subroutine adaptive_relax

!=======================================================================
!
!     AGILE: interpolation into changed independent variable
!
!=======================================================================

        subroutine adaptive_linint(new,status)

        integer, intent(out) :: status
        type(state_vector), intent(inout) :: new

!-----------------------------------------------------------------------
!
!     input:
!     adaptive_linint_template    state vector, no staggering assumed
!     new%y(new%jy(1:new%nq,1))   new grid point allocation
!
!     output:
!     new   state vector with all unknowns in old interpolated to new
!           grid point locations and this appropriately staggered
!     status   is set to 1 in case of a problem
!
!-----------------------------------------------------------------------

        integer :: nq,ny,nx,ixl,ixm,ixu,iax,ibx,iy,staggering
        integer, dimension(nnq) :: jq
        real, dimension(nnq) :: axt,ax,bxt,bx,f,y

!.....read from input file..............................................
        status = 0
        if (input_state%initial) call input_read_state
        ny = input_state%ny
        if (ny.ne.adaptive_linint_template%ny) then
          write(6,*) 'Error: adaptive_linint_template not properly set!'
          stop
        endif
        if (any(input_state%staggering(1:ny).eq.-1)) then
          staggering = -1
        elseif (any(input_state%staggering(1:ny).eq.1)) then
          staggering = 1
        else
          staggering = 0
        endif

!.....abbreviate notation and staggering................................
        nx = adaptive_linint_template%nq
        nq = new%nq
        call state_sumdx(adaptive_linint_template,axt)
        call state_sumdx(new,ax)
        if (staggering.eq.-1) then
          bxt(1) = axt(1)
          bxt(2:nq-1) = 0.5*(axt(2:nq-1)+axt(1:nq-2))
          bxt(nq) = axt(nq)
          bx(1) = ax(1)
          bx(2:nq-1) = 0.5*(ax(2:nq-1)+ax(1:nq-2))
          bx(nq) = ax(nq)
        elseif (staggering.eq.1) then
          bxt(1) = axt(1)
          bxt(2:nq-1) = 0.5*(axt(3:nq)+axt(2:nq-1))
          bxt(nq) = axt(nq)
          bx(1) = ax(1)
          bx(2:nq-1) = 0.5*(ax(3:nq)+ax(2:nq-1))
          bx(nq) = ax(nq)
        endif

!-----a-grid------------------------------------------------------------

!.....determine positions of new grid points............................
        do iax=1,nq
          ixl = 1
          ixu = nx
          do while (ixu-ixl.gt.1)
            ixm = (ixu+ixl)/2
            if ((axt(nx).gt.axt(1)).eqv.(ax(iax).gt.axt(ixm))) then
              ixl = ixm
            else
              ixu = ixm
            endif
          enddo
          jq(iax) = ixl
          f(iax) = (ax(iax)-axt(ixl))/(axt(ixl+1)-axt(ixl))
        enddo

!.....linearly interpolate variables....................................
        do iy=2,new%ny
          if (input_state%staggering(iy).eq.0) then
            y(1:nx) = adaptive_linint_template%y(
     1        adaptive_linint_template%jy(1:nx,iy))
            new%y(new%jy(1:nq,iy)) = y(jq(1:nq))
     1        + f(1:nq)*(y(jq(1:nq)+1)-y(jq(1:nq)))
          endif
        enddo

!-----b-grid------------------------------------------------------------

!.....determine positions of new grid points............................
        do ibx=1,nq
          ixl = 1
          ixu = nx
          do while (ixu-ixl.gt.1)
            ixm = (ixu+ixl)/2
            if ((bxt(nx).gt.bxt(1)).eqv.(bx(ibx).gt.bxt(ixm))) then
              ixl = ixm
            else
              ixu = ixm
            endif
          enddo
          jq(ibx) = ixl
          f(ibx) = (bx(ibx)-bxt(ixl))/(bxt(ixl+1)-bxt(ixl))
        enddo

!.....linearly interpolate variables....................................
        do iy=2,new%ny
          if (input_state%staggering(iy).ne.0) then
            y(1:nx) = adaptive_linint_template%y(
     1        adaptive_linint_template%jy(1:nx,iy))
            new%y(new%jy(1:nq,iy)) = y(jq(1:nq))
     1        + f(1:nq)*(y(jq(1:nq)+1)-y(jq(1:nq)))
          endif
        enddo

        end subroutine adaptive_linint

!=======================================================================
!
!     AGILE: adaptive grid equation
!
!=======================================================================

        subroutine adaptive_grid(jacobi,old,new,f,status)

        logical, intent(in) :: jacobi
        integer, intent(out) :: status
        real, dimension(nnc), intent(out) :: f
        type(state_vector), intent(in) :: old,new

!-----------------------------------------------------------------------
!
!     input:
!     jacobi   is true if the routine is called for derivatives
!     old   old state vector
!     new   new state vector
!
!     output:
!     f     residuum of grid equations,
!           filled in at positions (1+(iq-1)*ny, iq=2,nq-1)
!     status is set to 1 in case of a problem
!
!-----------------------------------------------------------------------
       
        logical, parameter :: smooth=.true.
        integer :: nq,nq1,i
        integer, dimension(nnq) :: jf
        real, dimension(2) :: bpsi
        real, dimension(nnq) :: psi

!.......initialisation..................................................
        status = 0
        nq = new%nq
        nq1 = nq-1
        if (input_allocation%initial) call input_read_allocation
        jf(2:nq1) = 1 + (/ (i, i=1,nq-2) /)*new%ny

!-----adaptive grid switched on-----------------------------------------
        if (input_allocation%grid) then
          call adaptive_deltaq(smooth,old,new,psi,bpsi,status)
          if (status.ne.0) return

!.....adaptive grid equations...........................................
          f(jf(2)) = bpsi(1)
          f(jf(3:nq-2)) = psi(3:nq-2) - psi(2:nq-3)
          f(jf(nq1)) = bpsi(2)

!-----adaptive grid switched off----------------------------------------
        else

!.....normal zones......................................................
          f(jf(2:nq-2)) = new%y(new%jy(2:nq-2,1))
     &      - old%y(old%jy(2:nq-2,1))

!.....allow equations_step to set zone nq1 OR n1........................
          f(jf(nq1)) = new%y(new%jy(nq,1))-new%y(new%jy(nq1,1))
        endif

        end subroutine adaptive_grid

!=======================================================================
!
!     AGILE: adaptive_deltaq
!
!=======================================================================

        subroutine adaptive_deltaq(smooth,old,new,psi,bpsi,status)

        logical, intent(in) :: smooth
        integer, intent(out) :: status
        real, dimension(2), intent(out) :: bpsi
        real, dimension(nnq), intent(out) :: psi
        type(state_vector), intent(in) :: old,new

!-----------------------------------------------------------------------
!
!     input:
!     smooth  if true, apply space smoothing
!     old     old state vector
!     new     new state vector
!
!     output:
!     psi   the grid is adapted if psi=delta(q) is constant for q=1:nq-1
!     bpsi  equations for boundary zones
!     status   is set to true if there is a problem with new state
!
!     we don't correct for the staggering in psi(1) and psi(nq-1)
!     because these values are only used in the initial allocation where
!     we always assume states without staggering.
!
!-----------------------------------------------------------------------

        logical :: space,time
        integer :: nq,nq1,nq2,ny
        integer :: iq,iy,staggering
        real :: f,taudt
        real, dimension(nny) :: a
        real, dimension(nnq) :: y,R,n,n0,n_tilde,n0_tilde,n_hat,tmp
        real, dimension(nnq) :: scl,relscl,dy

!.....initialise........................................................
        status = 0
        nq = new%nq
        nq1 = nq-1
        nq2 = nq-2
        ny = new%ny
        if (input_action%initial) call input_read_action
        if (input_allocation%initial) call input_read_allocation
        a(1:ny) = input_allocation%yloc(1:ny)
        f = input_allocation%alpha * (input_allocation%alpha + 1.)
        if (input_state%initial) call input_read_state
        if (any(input_state%staggering(1:ny).eq.-1)) then
          staggering = -1
        elseif (any(input_state%staggering(1:ny).eq.1)) then
          staggering = 1
        else
          staggering = 0
        endif

!-----calculating inverse generalized arc length------------------------

!.....first independent variable........................................
        iy=1
          
!.....scaling...........................................................
        call state_sumdx(new,y)
        if (input_allocation%scaling(iy).eq.2) then
          relscl(1:nq) = 1./sqrt(y(1:nq)**2 + old%ysclmin(iy)**2)
          scl(1:nq1) = 0.5*(relscl(2:nq) + relscl(1:nq1))
        else
          scl(1:nq1) = 1./old%yscl(iy)
        endif

!.....specific arc length...............................................
        y(1:nq) = new%y(new%jy(1:nq,iy))
        dy(1:nq1) = new%dx(2:nq) + y(2:nq) - y(1:nq1)
        R(1:nq1) = (a(iy)*scl(1:nq1)*dy(1:nq1))**2

!.....then loop over other variables....................................
        do iy=2,ny
          if (a(iy).gt.0.) then
          
!.....scaling...........................................................
            y(2:nq1) = new%y(new%jy(2:nq1,iy))
            if (input_allocation%scaling(iy).eq.2) then
              relscl(2:nq1) = 1./sqrt(y(2:nq1)**2 + old%ysclmin(iy)**2)
              scl(2:nq2) = 0.5*(relscl(3:nq1) + relscl(2:nq2))
            else
              scl(2:nq2) = 1./old%yscl(iy)
            endif

!.....specific arc length...............................................
            dy(2:nq2) = a(iy)*scl(2:nq2)*(y(3:nq1) - y(2:nq2))

!.....sum and invert....................................................
            R(2:nq2) = R(2:nq2) + dy(2:nq2)**2
          endif
        enddo
        R(1:nq1) = 1./sqrt(R(1:nq1))

!-----smoothing---------------------------------------------------------
        space = ((input_allocation%alpha.gt.0.).and.smooth)
        time = ((input_allocation%tau.gt.0.).and.(new%t.ne.old%t))

!.....new grid point concentration......................................
        y(1:nq) = new%y(new%jy(1:nq,1))
        tmp(1:nq1) = new%dx(2:nq) + y(2:nq) - y(1:nq1)
        if (any(tmp(1:nq1).eq.0.)) then
          if (input_action%message.ge.2)
     1      write(6,*) 'Warning: delta(x)=0 in adaptive_deltaq()!'
          status = 1
          return
        endif
        n(1:nq1) = 1./tmp(1:nq1)
        if (staggering.eq.-1) then
          n(nq1) = 0.5*n(nq1)
        elseif (staggering.eq.1) then
          n(1) = 0.5*n(1)
        endif

!.....old grid point concentration......................................
        if (time) then
          y(1:nq) = old%y(old%jy(1:nq,1))
          tmp(1:nq1) = old%dx(2:nq) + y(2:nq) - y(1:nq1)
          n0(1:nq1) = 1./tmp(1:nq1)
          if (staggering.eq.-1) then
            n0(nq1) = 0.5*n0(nq1)
          elseif (staggering.eq.1) then
            n0(1) = 0.5*n0(1)
          endif
        endif

!.....new space-smoothed grid point concentrations......................
        if (space) then
          n_tilde(1) = 1.
          n_tilde(2:nq2) = 1.
     1      - f*(n(3:nq1)/n(2:nq2) + n(1:nq2-1)/n(2:nq2) - 2.)
          n_tilde(nq1) = 1.

!.....old space-smoothed grid point concentrations......................
          if (time) then
            n0_tilde(1) = n0(1)/n(1)
            n0_tilde(2:nq2) = (n0(2:nq2)
     1        - f*(n0(3:nq1) + n0(1:nq2-1) - 2.*n0(2:nq2)))/n(2:nq2)
            n0_tilde(nq1) = n0(nq1)/n(nq1)
          endif

!.....not space-smoothed grid point concentrations......................
        else
          n_tilde(1:nq1) = 1.
          if (time) n0_tilde(1:nq1) = n0(1:nq1)/n(1:nq1)
        endif

!.....time-smoothed grid point concentration............................
        if (time) then
          taudt = input_allocation%tau/(new%t-old%t)
          n_hat(1:nq1) = n_tilde(1:nq1)
     1      + taudt*(n_tilde(1:nq1) - n0_tilde(1:nq1))
          bpsi(1) = n(2)-n(1)+taudt*(n(2)-n(1)-n0(2)+n0(1))
          bpsi(2) = n(nq1)-n(nq2)+taudt*(n(nq1)-n(nq2)-n0(nq1)+n0(nq2))

!.....not time-smoothed grid point concentration........................
        else
          n_hat(1:nq1) = n_tilde(1:nq1)
          bpsi(1) = n(2) - n(1)
          bpsi(2) = n(nq1) - n(nq2)
        endif
            
!.....apply smoothing...................................................
        psi(1:nq1) = R(1:nq1)*n_hat(1:nq1)

!.....equations for boundaries..........................................
          
        end subroutine adaptive_deltaq

!=======================================================================

      end module adaptive_grid_module

!***********************************************************************
