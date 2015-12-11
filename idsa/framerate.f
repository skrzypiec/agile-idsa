!=======================================================================
!
!     IDSA: framerate, set a frame around s and ye and evaluate rates
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine framerate(new,dt,d,slog,y,dr,ft,fs,etai,alpha,
     &  sg,yg,ug,dfdtg,slg,status)

      use egroup_module
      use species_module

      implicit none

      logical :: new
      integer :: status
      real, intent(in) :: dt,d,slog,dr
      real, dimension(nuc+1), intent(in) :: y
      real, dimension(ne,nut), intent(in) :: ft,fs,etai,alpha
      real, dimension(2) :: sg,yg
      real, dimension(2,2), intent(out) :: ug
      real, dimension(ne,nut,2,2), intent(out) :: dfdtg,slg

!-----------------------------------------------------------------------
!
!     ML, CITA (2004)
!     new ... flag to force update or end Newton-Raphson
!     dt ... time step [s]
!     d ... density [g/cm^3]
!     slog ... log entropy [kB/baryon]
!     y ... abundances [particles/baryon]
!     dr ... zone width
!     ft ... neutrino distribution functions of trapped neutrinos
!          ft = f(E) * 4.*pi/(h*c)**3*E**2*dE*mb/d  ![particles/baryon]
!     fs ... neutrino distribution functions of streaming neutrinos
!          fs = f(E) * 4.*pi/(h*c)**3*E**2*dE*mb/d  ![particles/baryon]
!     etai ... main factor in alpha for implicit finite differencing [1/s]
!     alpha ... diffusive neutrino source [particles/baryon/s]
!     sg ... log temperature on frame corners [kB/baryon]
!     yg ... electron fraction on frame corners [electrons/baryon]
!     ug ... specific internal energy on frame corners [MeV/baryon]
!     dfdtg ... f changes on frame corners [particles/baryon/s]
!     slg ... neutrino sources on frame corners [particles/baryon/s]
!
!-----------------------------------------------------------------------

      real, parameter :: fract=0.01,fracy=0.01

      integer :: i,j
      real :: seos,teos
      real, dimension(nuc+1) :: yeos,mu

!.....check input.......................................................
      status = 0.
      if (slog.gt.100.) then
        write(6,*) 'Warning: large slog in framerate.f'
        status = 21
        return
      endif

!.....check if a new frame is required..................................
      if (new.or.(sg(1).gt.slog).or.(slog.gt.sg(2)).or.
     &    (yg(1).gt.y(je)).or.(y(je).gt.yg(2))) then
        new = .true.
        seos = 10.**slog
        sg(1) = seos - fract*seos
        sg(2) = seos + fract*seos
        yg(1) = y(je) - fracy*y(je)
        yg(2) = y(je) + fracy*y(je)
     
!.....build new frame...................................................
        yeos = y
        do j=1,2
          yeos(je) = yg(j)
          do i=1,2
            seos = sg(i)

!.....define thermodynamic conditions...................................
            call eosinterf(d,seos,yeos,mu,teos,ug(i,j),status)
            if (status.ne.0) then
              write(6,*) 'Error: eos status=',status,' in framerate.f'
              write(6,*) 'status=1 --> density exceeds tabulated range'
              write(6,*) 'status=2 --> Ye exceeds tabulated range'
              write(6,*) 'status=3 --> entropy exceeds tabulated range'
              stop
            endif

!.....calculate rates...................................................
            call uydot(dt,d,teos,yeos,mu,dr,ft,fs,etai,alpha,
     &        dfdtg(1,1,i,j),slg(1,1,i,j))
          enddo
        enddo
        sg = log10(sg)
      endif

      end subroutine framerate

!=======================================================================
