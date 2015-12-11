!=======================================================================
!
!     IDSA: uyupdate, implicitly updates energy and electron fraction
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!     The Isotropic Diffusion Source Approximation was developed at the
!     University of Basel, funded by the Swiss National Science Foundation
!     under grant No. PP002-106627/1. It is documented in Liebendoerfer,
!     Whitehouse and Fischer, Astrophysical Journal 698 (2009) p. 1174.
!
!   Agile-IDSA is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License version 3
!   as published by the Free Software Foundation.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!=======================================================================

      subroutine uyupdate(dt,d,s,ye,dr,ft,fs,etai,alpha,iray,
     &  sout,yout,dfdt,sl,status)

      use egroup_module
      use frameint_module
      use species_module

      implicit none

      integer, intent(in) :: iray
      integer :: status
      real, intent(in) :: dt,d,s,ye,dr
      real, dimension(ne,nut), intent(in) :: ft,fs,etai,alpha
      real, intent(out) :: sout,yout
      real, dimension(ne,nut), intent(out) :: dfdt,sl

!-----------------------------------------------------------------------
!
!     Input:
!     dt ... time step [s]
!     d ... density [g/cm^3]
!     s ... entropy [kB/baryon]
!     ye ... electron fraction []
!     dr ... propagation length [cm]
!     ft ... neutrino distribution functions of trapped neutrinos
!          ft = f(E) * 4.*pi/(h*c)**3*E**2*dE*mb/d  ![particles/baryon]
!     fs ... neutrino distribution functions of streaming neutrinos
!          fs = f(E) * 4.*pi/(h*c)**3*E**2*dE*mb/d  ![particles/baryon]
!     etai ... main factor in alpha for implicit finite differencing [1/s]
!     alpha ... diffusive neutrino source [particles/baryon/s]
!     iray ... dimension array index (1 for 1D, 2 for 3D)
!
!     Output:
!     sout ... entropy after time step [kB/baryon]
!     yout ... electron fraction after time step [particles/baryon]
!     dfdt ... neutrino src from equil. adjustments [particles/baryon/s]
!     sl ... total approximate neutrino source [particles/baryon/s]
!
!-----------------------------------------------------------------------

      real, parameter :: eps=1.e-10

      logical :: new
      integer, parameter :: itmax=50
      integer :: it
      real, dimension(nuc+1) :: yeos
      real, dimension(ne,nut) :: flxe,dum3,dum4
      real :: u0,y0,det,dum1,dum2
      real, dimension(2) :: x,sg,yg,res,dx
      real, dimension(2,2) :: ug,jac
      real, dimension(ne,nut,2,2) :: dfdtg,slg

!.....construct energy groups...........................................
      if (egroup_empty) call egroup

!.....tabulate rates....................................................
      yeos(je) = ye
      x = (/log10(s),ye/)
      new = .true.
      call framerate(new,dt,d,x(1),yeos,dr,ft,fs,etai,alpha,
     &  sg,yg,ug,dfdtg,slg,status)
      if (status.ne.0) goto 100

!.....initial values....................................................
      call frameint(x(1),x(2),sg,yg,ug,u0,dum1,dum2)
      y0 = ye

!.....do Newton-Raphson iterations......................................
      it = 0
      do while (new)
        it = it+1
        if (it.gt.itmax) then
          status = 22
          goto 100
        endif

!.....calculate residuum and jacobian...................................
        call residuum(dt,u0,y0,sg,yg,ug,slg,x,res,jac)

!.....apply corrections.................................................
        det = jac(1,1)*jac(2,2) - jac(2,1)*jac(1,2)
        if (det.eq.0.) then
          write(6,*) 'Warning: det=0 in uyupdate.f'
	  write(6,*) x(1),x(2),sg,yg,ug
	  write(6,*) dt,u0,y0
	  write(6,*) slg,x,res
	  write(6,*) jac
          status = 23
          goto 100
        endif
        dx(1) = -(jac(2,2)*res(1) - jac(1,2)*res(2))/det
        dx(2) = -(jac(1,1)*res(2) - jac(2,1)*res(1))/det
        x = x+dx
        
!.....check if new tabulation is required and if converged..............
        yeos(je) = x(2)
        new = .false.
        call framerate(new,dt,d,x(1),yeos,dr,ft,fs,etai,alpha,
     &    sg,yg,ug,dfdtg,slg,status)
        if (status.ne.0) goto 100
        if (.not.new) new = (abs(dx(1))+abs(dx(2)/x(2)).gt.eps)
      enddo

!.....set output entropy and electron fraction..........................
      sout = 10.**x(1)
      yout = x(2)

!.....interpolate dfdt and sl.........................................
      call frameint(x(1),x(2),sg,yg,dfdtg,dfdt,dum3,dum4)
      call frameint(x(1),x(2),sg,yg,slg,sl,dum3,dum4)

!.....regular return....................................................
      return

!.....return on error...................................................
100   continue
      sout = s
      yout = ye
      dfdt = 0.
      sl = 0.

      end subroutine uyupdate

!=======================================================================
