!***********************************************************************
!
!     IDSA: spectrum_module
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module spectrum_module

      use egroup_module
      use species_module

      implicit none

      logical, parameter :: spectrum_disable=.true.
      integer, parameter :: nray=1 !number of rays
      integer :: nzone    !maximum number of zones

!.....energy-dependent corrections to Fermi-Dirac spectrum..............
      real, allocatable, dimension(:,:,:,:) :: spectrum_df
      real, allocatable, dimension(:,:,:,:) :: src ![particles/baryon/s]
      
      contains

!=======================================================================
!
!     IDSA: spectrum_initialise
!
!=======================================================================

        subroutine spectrum_initialise
	
        use input_file_module

        if (input_state%initial) call input_read_state
        nzone = input_state%nq
	allocate(spectrum_df(ne,nut,nzone,nray))
        allocate(src(ne,nut,nzone,nray))
	spectrum_df = 0.
        src = 0.
	   
	end subroutine spectrum_initialise

!=======================================================================
!
!     IDSA: spectrum_update
!
!=======================================================================

        subroutine spectrum_update(i,dt,d,ynu,znu,dfdt,iray)

        use units_module

        implicit none

        integer, intent(in) :: i,iray
        real, intent(in) :: dt,d
        real, dimension(nut), intent(in) :: ynu,znu
        real, dimension(ne,nut), intent(in) :: dfdt

!-----------------------------------------------------------------------
!
!     Input:
!     i ... zone identifier
!     dt ... time step [s]
!     d ... density [g/cm^3]
!     ynu ... neutrino fraction [particles/baryon]
!     znu ... neutrino specific energy [MeV/baryon]
!     dfdt ... neutrino src from equil. adjustments [particles/baryon/s]
!     iray ... array index for dimension, 1 for 1D, 2 for 3D
!
!     Output:
!     spectrum_df ... residual spectrum []
!
!-----------------------------------------------------------------------

        integer :: ie,inu,status
        real, parameter :: fac=4.*units%pi/(units%h*units%c)**3
        real :: nb,nd,ed,beta,eta,fermi,ynutmp,znutmp
        real :: a,b,ae,be,be2
        real, dimension(ne) :: ps,f0,f,df

!.......return immediately if spectral corrections are disabled.........
        if (spectrum_disable) then
          spectrum_df = 0.
          return
        endif

!.....initialize........................................................
        nb = d/units%mb           !baryon density
        ps = fac*E**2*dE*units%mb !particles/cm^3 * g/baryon
        if (i.gt.nzone) then
          write(6,*) 'Error: i>nzone in spectrum_module.f'
	  write(6,*) 'i: ',i,' nzone:',nzone
          stop
        endif
        status = 0

!.....loop over neutrino species........................................
        do inu=1,nut

!.....old distribution function.........................................
          nd = nb*ynu(inu)        !neutrino number density
          ed = nb*znu(inu)        !neutrino specific energy
          call etafit(nd,ed,beta,eta,status)
          if (status.ne.0) then
            spectrum_df(:,inu,i,iray) = 0.
            cycle
          endif
          do ie=1,ne
            f0(ie) = fermi(beta*E(ie)-eta)
          enddo

!.....new ynu,znu.......................................................
          ynutmp = ynu(inu) + sum(dfdt(:,inu))*dt   !particles/baryon
          znutmp = znu(inu) + sum(dfdt(:,inu)*E)*dt !MeV/baryon

!.....new distribution function.........................................
          nd = nb*ynutmp        !neutrino number density
          ed = nb*znutmp        !neutrino specific energy
          call etafit(nd,ed,beta,eta,status)
          if (status.ne.0) then
            spectrum_df(:,inu,i,iray) = 0.
            cycle
          endif
          do ie=1,ne
            f(ie) = fermi(beta*E(ie)-eta)
          enddo

!.....spectral changes..................................................
          df = f0 + spectrum_df(:,inu,i,iray) + dfdt(:,inu)*dt*d/ps - f

!.....enforce zero net number and energy of spectral corrections........
          a = sum(df*ps)
          b = sum(f*ps)
          ae = sum(df*E*ps)
          be = sum(f*E*ps)
          be2 = sum(f*E*E*ps)
          spectrum_df(:,inu,i,iray) = df
     &      + (a*be2-ae*be - (a*be-ae*b)*E)/(be*be - b*be2)*f
        enddo

        end subroutine spectrum_update

!=======================================================================
!
!     IDSA: spectrum_add
!
!=======================================================================

        subroutine spectrum_add(i,d,ynu,znu,iray,f,status)

        use units_module

        implicit none

        integer, intent(in) :: i,iray
        real, intent(in) :: d
        real, dimension(nut), intent(in) :: ynu,znu
        real, dimension(ne,nut), intent(out) :: f
        integer, dimension(nut), intent(out) :: status

!-----------------------------------------------------------------------
!
!     Input:
!     i ... zone identifier
!     d ... density [g/cm^3]
!     ynu ... neutrino fraction []
!     znu ... neutrino specific energy [MeV/baryon]
!     iray ... array index for dimensions, 1 for 1D, 2 for 3D
!
!     Output:
!     f ... distribution function with corrected spectrum []
!     status ... 0 if no problems occured
!
!-----------------------------------------------------------------------

        integer :: ie,inu
        real :: nb,nd,ed,beta,eta,fermi
        real, dimension(ne) :: f0

!.....initialize........................................................
        nb = d/units%mb           !baryon density
        if (i.gt.nzone) then
          write(6,*) 'Error: i>nzone in spectrum_module.f add'
	  write(6,*) 'i: ',i,' nzone:',nzone
          stop
        endif

        status = 0
       
!.....loop over neutrino species........................................
        do inu=1,nut

!.....old distribution function.........................................
          nd = nb*ynu(inu)        !neutrino number density
          ed = nb*znu(inu)        !neutrino specific energy
          call etafit(nd,ed,beta,eta,status(inu))
          if (status(inu).ne.0) then
            f0 = 0.
          else
            do ie=1,ne
              f0(ie) = fermi(beta*E(ie)-eta)
            enddo
          endif

!.....add spectral corrections to distribution function.................
          f(:,inu) = f0 + spectrum_df(:,inu,i,iray)
        enddo

        end subroutine spectrum_add

!=======================================================================

      end module spectrum_module

!***********************************************************************
