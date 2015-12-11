!***********************************************************************
!
!     AGILE: newton_module
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module newton_module

      use adaptive_grid_module
      use agile_parameter_module
      use input_file_module
      use sparse_matrix_module
      use state_vector_module

      implicit none

!=======================================================================
!
!     AGILE: definition of correction control flags
!
!=======================================================================

      type newton_control_flags
        logical :: block       !setup of block-diagonal dummy jacobian
        logical :: analyze     !find nonzeros in jacobian
        logical :: jac_update  !calculate a new jacobian
        logical :: lu_update   !update band-diagonal LU-decomposition
        logical :: banddiag    !solve system by back-substitution
        logical :: bigrad      !solve system by bi-conjugate gradient
        real :: tol            !tolerance for linear system
      end type newton_control_flags
      type(newton_control_flags) :: newton_flag

!-----------------------------------------------------------------------
!
!     if one of these flags is raised then the corresponding action
!     is executed at the next occasion. The flag is set back to .false.
!
!-----------------------------------------------------------------------
        
!=======================================================================
!
!     AGILE: interface for physical equations
!
!=======================================================================

      interface
	subroutine equations(jacobi,old,new,f,status)
        use state_vector_module
        logical, intent(in) :: jacobi
        integer, intent(out) :: status
        real, dimension(nnc), intent(out) :: f
        type(state_vector), intent(in) :: old,new
        end subroutine equations
      end interface

      contains

!=======================================================================
!
!     AGILE: newton_raphson
!
!=======================================================================

	subroutine newton_raphson(itmax,slowc,conv,
     1    old,new,status,equations)

        external equations

        integer, intent(out) :: status
        integer, intent(in) :: itmax
        real, intent(in) :: slowc,conv
        type(state_vector), intent(in) :: old
        type(state_vector), intent(inout) :: new

!-----------------------------------------------------------------------
!
!     input:
!     itmax       maximal number of iterations
!     slowc       damping in case of overshooting
!     conv        defining convergence
!     old         old state
!     new         guess for solution
!     equations   equations to be solved
!     set the following newton_flags before calling this routine:
!     newton_flag%block
!     newton_flag%analyze
!     newton_flag%banddiag
!     newton_flag%bigrad
!
!     output:
!     new         improved guess
!     status      status=0 if no problem was encountered
!
!-----------------------------------------------------------------------

        integer :: nq,ny,n,it,iy,iq
        integer, dimension(1) :: iloc
        integer, dimension(nny) :: icor
        real :: err,damp
        real, dimension(1) :: rval
        real, dimension(nnq) :: aloc
        real, dimension(nnc) :: f,cor

!.....initialise........................................................
        status = 0
        nq = new%nq
        ny = new%ny
        n = nq*ny
        newton_flag%tol = conv/10.
        newton_flag%lu_update = .true.
        if (any(old%yscal(1:n).le.0.)) then
          write(6,*) 'Error: old state not scaled in newton_raphson()'
          stop
        endif

!.....adjust zone differences to new shift..............................
        call state_newdx(new)

!.....iteration loop....................................................
	do it=1,itmax

!.....calculate corrections.............................................
          newton_flag%jac_update = .true.
          if (newton_flag%banddiag) newton_flag%lu_update = .true.
          call newton_correction(old,new,f,cor,status,equations)
          if (status.ne.0) then
            write(6,*) 'Warning: Problem occured in newton_correction'
            goto 100
          endif
          rval = maxval(abs(cor(1:n)),1)
          err = rval(1)

!.....write errors......................................................
          if (input_action%message.ge.2) then
            write(6,11)
            write(6,22) 'variable','grid point','correction','scale'
            do iy=1,ny
              aloc(1:nq) = abs(cor(new%jy(1:nq,iy)))
              iloc = maxloc(aloc(1:nq))
              iq = iloc(1)
              write(6,33) iy,iq,
     1          cor(new%jy(iq,iy)),old%yscal(old%jy(iq,iy))
            enddo
            write(6,44) it,err
          endif
11    format(72('-'))
22    format(4a12)
33    format(2i12,2g12.4)
44    format(1x,'Newton-Raphson iteration: ',i4,', residuum: ',g12.4)

!.....investigating errors..............................................
          if (err.gt.10.) then
            if (input_action%message.ge.1)
     1        write(6,*) 'Info: guess rejected by newton_raphson()'
            status = 1
            goto 100
          endif

!.....pause.............................................................
          if (input_action%message.ge.4) then
            write(6,*) 'Info: this pause allows you to check the file'
            write(6,*) 'jacobian.tmp (press return to continue)'
            read(5,*)
          endif

!.....apply corrections.................................................
          damp = slowc/max(slowc,err)
          new%y(1:n) = new%y(1:n) + damp*cor(1:n)*old%yscal(1:n)

!.....if converged return with solution.................................
          if (err.lt.conv) return
        enddo   !iteration loop

!.....cleanup the guess and return without solution.....................
        if (input_action%message.ge.1)
     1    write(6,*) 'Warning: no convergence in newton_raphson()'
100     continue
        new%y(new%jy(1:nq,1)) = 0.
        status = 1

        end subroutine newton_raphson

!=======================================================================
!
!     AGILE: newton_correction
!
!=======================================================================

        subroutine newton_correction(old,new,f,cor,status,equations)

        external equations

        integer, intent(out) :: status
        real, dimension(nnc), intent(out) :: f,cor
        type(state_vector), intent(in) :: old
        type(state_vector), intent(inout) :: new

!-----------------------------------------------------------------------
!
!       input:
!       old   last time state vector or dummy in static case
!       new   guess state vector
!       equations   subroutine with the equations to solve
!       all newton_flags have to be set before calling this routine!
!
!       output:
!       f   residua
!       cor   (one iteration)-corrections to guess
!       status   set to 1 on output if new is to bad
!
!-----------------------------------------------------------------------

        logical, parameter :: jacobi=.false.
        integer :: i,iter
        integer, save :: m1,m2
        real :: err
        type(sparse_matrix), save :: a

!-----calculate residua-------------------------------------------------

!.....evaluate physical equations with new............................
        status = 0
        if (input_action%message.ge.3) write(6,11)
        call state_fill_zones(new,status)
        if (status.ne.0) return
	call equations(jacobi,old,new,f,status)
	if (status.ne.0) return

!.....evaluate adaptive grid equation...................................
        if (input_action%message.ge.3) write(6,22)
        call adaptive_grid(jacobi,old,new,f,status)
        if (status.ne.0) return

!-----calculate jacobian------------------------------------------------
        if (newton_flag%jac_update) then

!.....initial setup of jacobian.........................................
          if (newton_flag%block) then
            if (input_jacobi%initial) call input_read_jacobi
            if (input_action%message.ge.3) write(6,33)
            call sparse_block(new%nq,new%ny,input_jacobi%np,
     1        input_jacobi%voffset,input_jacobi%hoffset,a)
            if (input_action%message.ge.3) write(6,44)
            call sparse_analyze(a,m1,m2)
          endif

!.....calculate jacobian................................................
          if (input_action%message.ge.3) write(6,55)
	  call newton_jacobian(old,new,f,a,status,equations)
          if (status.ne.0) return
          if (input_action%message.ge.3)
     1      call sparse_write('sparse.tmp',a)

!.....eliminate nonzeros in jacobian, but not the first time............
          if (newton_flag%analyze) then
            if (newton_flag%block) then
              newton_flag%block = .false.
            else
              if (input_action%message.ge.3) write(6,66)
              call sparse_analyze(a,m1,m2)
              newton_flag%analyze = .false.
            endif
          endif
          newton_flag%jac_update = .false.
        endif

!-----LU-decomposition--------------------------------------------------
        if (input_action%message.ge.3) write(6,77)
	call sparse_solve(a,f,cor,status,input_action%message,m1,m2)
	if (status.ne.0) return

11      format(1x,'Info: evaluating equations')
22      format(1x,'Info: evaluating adaptive grid')
33      format(1x,'Info: building initial jacobian')
44      format(1x,'Info: analyzing initial jacobian')
55      format(1x,'Info: calculating jacobian')
66      format(1x,'Info: analyzing jacobian')
77      format(1x,'Info: decomposing band-diagonal system')
88      format(1x,'Info: preconditioning sparse system')
99      format(1x,'Info: iterations, error ',i12,g12.4)
111     format(1x,'Error: no solver specified in newton_correction.f!')
122     format(1x,'Info: solving band-diagonal system')
133     format(1x,'Info: solving sparse system')

        end subroutine newton_correction

!=======================================================================
!
!     AGILE: newton_jacobian
!
!=======================================================================

	subroutine newton_jacobian(old,new,f,a,status,equations)

        external equations

        integer, intent(out) :: status
	type(state_vector), intent(in) :: old,new
	type(sparse_matrix), intent(inout) :: a
	real, dimension(nnc), intent(in) :: f

!-----------------------------------------------------------------------
!
!     input:
!     old      state vector at time t
!     new      state vector at time t+dt
!     f        residua of equations evaluated with oldstate and new
!     a        sparse matrix describing nonzero positions for jacobian
!              and providing coupling information
!
!     output:
!     a        jacobian as sparse matrix
!     status   is set to 1 if variation in guess is bad
!
!     remark:
!     The jacobian is calculated by numerical derivatives. The coupling
!     information of the sparse matrix allows the separation of this
!     process into independent groups. The different groups do not
!     share data, so that they may be calculated in parallel. Each
!     group itself involves only evaluations of the equations over the
!     whole computational domain, so that equations() can be vectorized
!     for each single group evaluation. The parallelisation is not
!     implemented (yet?).
!
!-----------------------------------------------------------------------

        logical, parameter :: jacobi=.true.
	integer :: ny,nq,n,ig,ib,ib0,ib1,ia,ia0,ia1
	real, parameter :: eps=1.e-8
	real, dimension(nnc) :: tmpf,h,yvar
	type(state_vector) :: tmp

        status = 0
	nq = new%nq
	ny = new%ny
	n = nq*ny
	
!.....prepare variations................................................
        yvar(1:n) = new%y(1:n) + eps*old%yscal(1:n)
        h(1:n) = (yvar(1:n) - new%y(1:n))/old%yscal(1:n)
        if (any((h(1:n).eq.0.))) then
          write(6,*) 'h equal 0 in newton_jacobi()!'
          stop
        endif

!.....loop over groups..................................................
	do ig=1,a%ng

!.....make a copy of new...........................................
	  tmp = new
	  
!.....select columns belonging to actual group..........................
	  ib0 = a%jb(ig)+1
	  ib1 = a%jb(ig+1)

!.....apply variation to variables that correspond to columns jc........
          tmp%y(a%jc(ib0:ib1)) = yvar(a%jc(ib0:ib1))

!.....evaluate physical equations.......................................
          call state_fill_zones(tmp,status)
          if (status.ne.0) return
	  call equations(jacobi,old,tmp,tmpf,status)
          if (status.ne.0) return

!.....evaluate adaptive grid equation...................................
          call adaptive_grid(jacobi,old,tmp,tmpf,status)
          if (status.ne.0) return

!.....loop over columns.................................................
	  do ib=ib0,ib1

!.....calculate jacobian................................................
	    ia0 = a%ja(a%jc(ib))+1
	    ia1 = a%ja(a%jc(ib)+1)
            do ia=ia0,ia1
              if (a%numeric(ia)) a%a(ia) = (tmpf(a%jr(ia))
     1          - f(a%jr(ia)))/h(a%jc(ib))
            enddo
	  enddo   !loop over columns
	enddo   !loop over groups

	end subroutine newton_jacobian

!=======================================================================

      end module newton_module

!***********************************************************************
