!***********************************************************************
!
!     AGILE: sparse_matrix_module
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module sparse_matrix_module

      use agile_parameter_module
      use input_file_module

      implicit none

!=======================================================================
!
!     AGILE: type definition for sparse_matrix
!
!=======================================================================

      type sparse_matrix

!.....matrix definition.................................................
        integer :: nc                  !number of columns
                                       !begin row attributes a:
        real, dimension(nns) :: a      !nonzero coefficients of matrix
        logical, dimension(nns) :: numeric   !mode for derivatives
                                       !end row attributes.
        integer, dimension(nns) :: jr  !the row attribute ia belongs
                                       !to row jr(ia)
        integer, dimension(nnc+1) :: ja   !the attributes ia belonging
                                       !to column ic are at position
                                       !ja(ic) < ia <= ja(ic+1)

!.....a group contains columns with disjunkt coupling to rows...........
        integer :: ng                  !number of groups
                                       !begin column attributes b:
                                       !end column attributes.
        integer, dimension(nnc) :: jc  !the column attribute ib belongs
                                       !to column jc(ib)
        integer, dimension(nnc+1) :: jb   !the attributes ib belonging
                                       !to group ig are at position
                                       !jb(ig) < ib <= jb(ig+1)
      end type sparse_matrix

      contains

!=======================================================================
!
!     AGILE: create blockdiagonal sparse matrix
!
!=======================================================================

        subroutine sparse_block(nq,ny,ncouple,voffset,hoffset,a)

        integer :: nq,ny,ncouple,voffset,hoffset
        type(sparse_matrix) :: a

!-----------------------------------------------------------------------
!
!     input:
!     nq   number of grid points
!     ny   number of variables per grid point
!     ncouple   number of coupled zones
!     voffset   vertical offset of block-diagonal structure
!     hoffset   horizontal offset of block-diagonal struture
!
!     output:
!     a   sparse, block-diagonal matrix with all coefficients=1
!    
!-----------------------------------------------------------------------

        integer :: ic,ia,iq,iy,ir0,ir1,ir,m1,m2

!.....check validity of input...........................................
        a%nc = nq*ny
        if (a%nc.gt.nnc) then
          write(6,*) 'nq*ny larger than nnc in sparse_block()!'
          stop
        endif
        if ((voffset.lt.0).or.(hoffset.lt.0)) then
          write(6,*) 'invalid offset in sparse_block()!'
          stop
        endif
        if (voffset.ge.hoffset*ny) then
          write(6,*) 'voffset too large in sparse_block()!'
          stop
        endif
        if (hoffset.ge.ncouple) then
          write(6,*) 'hoffset larger than ncouple in sparse_block()!'
          stop
        endif
        
!.....initialise........................................................
        ic = 0
        ia = 0

!.....loop over grid points.............................................
        do iq=1,nq

!.....loop over variables...............................................
          do iy=1,ny

!.....create new column.................................................
            ic = ic+1
            a%ja(ic) = ia

!.....loop over rows....................................................
            ir0 = (iq+hoffset - ncouple)*ny - voffset + 1
            ir1 = (iq+hoffset)*ny - voffset
            do ir=max(ir0,1),min(ir1,a%nc)

!.....set matrix........................................................
              ia = ia+1
              if (ia.gt.nns) then
                write(6,*) 'ia exceeds nns in sparse_block()!'
                stop
              endif
              a%a(ia) = 1.
              a%numeric(ia) = .true.
              a%jr(ia) = ir
            enddo   !loop over rows
          enddo   !loop over variables
        enddo   !loop over grid points
        a%ja(ic+1) = ia

!.....set coupling information..........................................
        call sparse_analyze(a,m1,m2)

        end subroutine sparse_block

!=======================================================================
!
!     AGILE: analyze sparse matrix
!
!=======================================================================

        subroutine sparse_analyze(a,m1,m2)

        integer, intent(out) :: m1,m2
        type(sparse_matrix), intent(inout) :: a

!-----------------------------------------------------------------------
!
!     input:
!     a   sparse matrix, eventually containing zero coefficients and
!         an undefined coupling information
!
!     output:
!     a   same matrix, zeros eliminated, coupling information provided
!     m1   number of subdiagonals
!     m2   number of superdiagonals
!
!-----------------------------------------------------------------------

	logical :: finished,found
	logical, dimension(nnc) :: col_todo,rblocked,rexclude
	integer :: ic,ir,ia,iatmp,iastart,ib,ig

        if (input_action%initial) call input_read_action

!-----eliminate zeros---------------------------------------------------
        iatmp = 0
        m1 = 0
        m2 = 0

!.....loop over columns.................................................
        do ic=1,a%nc
          iastart = a%ja(ic)+1
          a%ja(ic) = iatmp

!.....loop over rows....................................................
          do ia=iastart,a%ja(ic+1)
            if (a%a(ia).ne.0.) then
              iatmp = iatmp+1
              a%a(iatmp) = a%a(ia)
              a%numeric(iatmp) = a%numeric(ia)
              ir = a%jr(ia)
              a%jr(iatmp) = ir
              if (ir.gt.ic) then
                m1 = max(m1,ir-ic)
              else
                m2 = max(m2,ic-ir)
              endif
            endif
          enddo   !loop over rows
        enddo   !loop over columns
        a%ja(a%nc+1) = iatmp

!-----update coupling information---------------------------------------

!.....initialise........................................................
	finished = .false.
        ig = 0
        ib = 0
	col_todo(1:a%nc) = .true.

!.....loop until all columns calculated.................................
	do while (.not.finished)
	  finished = .true.
          ig = ig+1
          a%jb(ig) = ib

!.....free all rows.....................................................
          rblocked(1:a%nc) = .false.
          rexclude(1:a%nc) = .false.

!.....loop over columns.................................................
	  do ic=1,a%nc
	    if (col_todo(ic)) then
	    
!.....check for forbidden rows coupled to that column...................
              do ia=a%ja(ic)+1,a%ja(ic+1)
                if (rexclude(a%jr(ia))) then
                  found = a%numeric(ia)
                else
                  found = rblocked(a%jr(ia))
                endif
                if (found) exit
              enddo

!.....forbidden row found: go to next column and try this one later.....
	      if (found) then
		finished = .false.

!.....no forbidden row found: add column to group and forbid rows.......
              else
                ib = ib+1
                a%jc(ib) = ic
		do ia=a%ja(ic)+1,a%ja(ic+1)
                  if (a%numeric(ia)) then
		    rblocked(a%jr(ia)) = .true.
                  else
		    rexclude(a%jr(ia)) = .true.
                  endif
		enddo
                col_todo(ic) = .false.
	      endif   !found
	    endif   !col_todo
	  enddo   !loop over columns
	enddo   !while not finished

!.....set endpoint of column list.......................................
        a%jb(ig+1) = ib
        a%ng = ig
        if (input_action%message.ge.3) write(6,11) a%ng,a%ja(a%nc+1)
11      format(1x,'Info: ',i6,' groups,',i12,' nonzeros')

	end subroutine sparse_analyze

!=======================================================================
!
!     AGILE: sparse_write
!
!=======================================================================

        subroutine sparse_write(filename,a)

        character(len=*) :: filename
        type(sparse_matrix) :: a

!-----------------------------------------------------------------------
!
!     input:
!     filename
!     a   sparse matrix
!
!-----------------------------------------------------------------------

        integer :: ig,ib,ic,ia

!.....open file.........................................................
        open(unit=1, file=filename)

!.....write coupling information........................................
        write(1,11) 'index','group','column'
        do ig=1,a%ng
          do ib=a%jb(ig)+1,a%jb(ig+1)
            write(1,22) ib,ig,a%jc(ib)
          enddo
        enddo

!.....write matrix......................................................
        write(1,11) 'index','column','row','coefficient','numeric'
        do ic=1,a%nc
          do ia=a%ja(ic)+1,a%ja(ic+1)
            write(1,33) ia,ic,a%jr(ia),a%a(ia),a%numeric(ia)
          enddo
        enddo

!.....close file and format.............................................
        close(unit=1)
11      format(3a12,a24,a12)
22      format(3i12)
33      format(3i12,g24.16,l12)

        end subroutine sparse_write
       
!=======================================================================
!
!     AGILE: sparse_solve
!
!=======================================================================

        subroutine sparse_solve(a,b,x,status,message,m1,m2)

        integer :: status
        integer, intent(in) :: message
        integer, intent(in) :: m1          !number of subdiagonal rows
        integer, intent(in) :: m2          !number of superdiagonal rows
        real, dimension(nnc), intent(in) :: b
        real, dimension(nnc), intent(out) :: x
        type(sparse_matrix) :: a

!-----------------------------------------------------------------------
!
!     input:
!     a    sparse matrix
!     b    right hand side of linear system
!     m1   number of subdiagonals
!     m2   number of superdiagonals
!
!     output
!     x = -(A^t)^(-1)*b
!     status= 0, exec. was successful
!     status=-i, ith parameter was illegal
!     status= i, u(i,i) was zero
!
!-----------------------------------------------------------------------

        integer :: nb,n                   !rank of matrix
        real, dimension(nnc) :: rmax      !row scaling for pivoting
        real, dimension(:,:), allocatable :: ab  !band-diagonal matrix
        integer, dimension(nnc) :: ipiv   !pivot index: row i was
                                          !interchanged with row ipiv(i)

        n = a%nc

!-----factorize band-diagonal approximation-----------------------------
        nb = 2*m1+m2+1
        allocate(ab(nb,n))
        call sparse_to_band(a,m1,m2,ab,nb,n,rmax)

!.....LU factorization of band matrix..................................
        call dgbtrf(n,n,m1,m2,ab,nb,ipiv,status)
        if (status.ne.0) return

!-----solve approximate system-----------------------------------------
        x(1:n) = -b(1:n)/rmax(1:n)
        
!.....direct system.....................................................
        call dgbtrs('N',n,m1,m2,1,ab,nb,ipiv,x,nnc,status)

        end subroutine sparse_solve

!=======================================================================
!
!     AGILE: convert sparse_matrix to band_matrix
!
!=======================================================================

        subroutine sparse_to_band(a,m1,m2,b,nb,n,rmax)

        integer, intent(in) :: m1,m2,nb,n
        type(sparse_matrix), intent(in) :: a
        real, dimension(nb,n) :: b  !band-diagonal matrix
        real, dimension(nnc) :: rmax !row scaling for pivoting

!-----------------------------------------------------------------------
!
!     input:
!     a   sparse matrix
!     m1   number of subdiagonals
!     m2   number of superdiagonals
!
!     output:
!     b   band-diagonal matrix, columns scaled
!           by maximal column element Cscl
!    
!-----------------------------------------------------------------------

        integer :: ic,is,s1,s2
        real :: tmp

!.....initialise band matrix............................................
        b = 0.
        rmax = 0.
     
!.....row maxima........................................................
        do ic=1,n
          s1 = a%ja(ic)+1
          s2 = a%ja(ic+1)
          do is = s1,s2
            tmp = abs(a%a(is))
            if (tmp.gt.rmax(a%jr(is))) rmax(a%jr(is)) = tmp
          enddo
        enddo

!.....loop over columns.................................................
        do ic = 1,n
          s1 = a%ja(ic)+1
          s2 = a%ja(ic+1)
          b(a%jr(s1:s2)+m2-ic+1+m1,ic) = a%a(s1:s2)/rmax(a%jr(s1:s2))
        enddo

        end subroutine sparse_to_band

!=======================================================================

      end module sparse_matrix_module

!***********************************************************************
