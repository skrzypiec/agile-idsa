!***********************************************************************
!
!     LEAK: leakage scheme for mu and tau neutrinos
!     Copyright: Albino Perego, 30.09.2011, GNU public license.
!
!***********************************************************************

      module mu_tau_leakage_module

      use egroup_module
      use units_module

      implicit none

      contains

!=======================================================================
!
!     LEAK: leakage scheme
!
!=======================================================================

      subroutine leakage(d,T,Ye,mfpm,opdepm,qloc,qdif,qeff)

!-----------------------------------------------------------------------
! This subroutine implements the leakage scheme for the energy
! release due to mu and tau neutrinos. It calculates,for every point 
! of the grid, the effective energy emission rate, expressed as a smooth 
! interpolation between the diffusion rate and the local production rate.
! For reference, see Rosswog&Liebendorfer(2003), Appendix
!-----------------------------------------------------------------------

      implicit none

      real, intent (in)  :: d      
      real, intent (in)  :: T      
      real, intent (in)  :: Ye
      real, dimension(ne), intent (in)  :: mfpm
      real, dimension(ne), intent (in)  :: opdepm

      real, intent (out) :: qloc   
      real, intent (out) :: qdif  
      real, intent (out) :: qeff 

!.....................................................................
!  Input and output parameters declaration                             
!                                                                     
!  The input terms are:                                               
!  d     ---> local density [g/cm^3]                                       
!  T     ---> local temperature [MeV]                                    
!  Ye    ---> local electron fraction
!  mfpm  ---> vector of local mu/tau neutrinos mean free path [cm]
!  opdepm --> vector of local mu/tau neutrinos optical depth
!                      
!  The output terms are
!  qloc ---> the local energy production rate [ergs/cm^3/s] 
!  qdif ---> the local energy diffusion rate [ergs/cm^3/s] 
!  qeff ---> the local effective energy emission rate [ergs/cm^3/s] 
!.....................................................................

      integer :: status
      real, dimension(ne) :: tdif
      real :: stat_weight,tmp

!.....check thermodynamic input
      if ((d.le.0.).or.(T.le.0.).or.(Ye.le.0.)) then
        write(6,*) 'Warning: Invalid input, mu_tau_leakage_module.f'
        write(6,*) d,T,Ye
        qloc = 0.
        qdif = 0.
        qeff = 0.
        return
      endif

!.....energy production rate
      call pair(d,t,Ye,qloc)

!.....diffusion timescale
      call dif_time(opdepm,mfpm,tdif)

!.....energy diffusion rate
      stat_weight = 4.
      call diffusion(t,tdif,stat_weight,qdif,status)

!.....interpolation between the production and the diffusion energy rate
      if (status.eq.0) then
        tmp = qloc + qdif
        if (tmp.eq.0.) then
          qeff = 0.
        else
          qeff = qloc*qdif/tmp
        endif
      else
        qeff = qloc
      endif
      if (isnan(qeff)) then
        write(6,*) "Error: qeff is NaN in mu_tau_leakage_module.f"
        write(6,*) qloc,qdif,tmp
        stop
      endif
 
      end subroutine leakage

!====================================================================
!
!     LEAK: Pair production rate
!
!====================================================================

      subroutine pair(rho,t,ye,qep)

!====================================================================
! This subroutine calculates the energy emission rate from the pair 
! production due to positron-electron annihilation. 
! To compute the rate of this process, we use the fitting formula 
! from Itoh et al.(1996), isolating the mu and tau contribution 
! from the total rate
!===================================================================

      implicit none

      real, intent(in) :: rho
      real, intent(in) :: t
      real, intent(in) :: ye

      real, intent(out) :: qep

!....................................................................
!     Input variables
!     rho ----> local density [g/cm^3]
!     t   ----> local temperature [MeV]
!     ye  ----> local electron fraction
!
!     Output variables
!     qep ----> energy production rate [ergs/g/s]
!....................................................................

!.....parameters and constants
      real, parameter :: cv = 0.5+2.*units%sinsqthetw ! vector coupling constant
      real, parameter :: ca = 0.5 ! axial coupling constant
      real, parameter :: pair_num = 1.

      real :: pref     ! dimensional factor
      real :: memw     ! mean molecular weight per electron (1/Ye)
      real :: lambda   ! dimensionless T parameter from Itoh et al. (1996)
      real :: xi       ! dimensionless parameter from Itoh et al. (1996)
      real :: a0,a1,a2 ! dimensionless coefficients from Itoh et al. (1996)
      real :: b1,b2,b3 ! dimensionless coefficients from Itoh et al. (1996)
      real :: cc       ! dimensionless coefficient from Itoh et al. (1996)
      real :: qpair    ! function from Itoh et al. (1996)
      real :: g        ! function from Itoh et al. (1996)
      real :: fpair    ! function from Itoh et al. (1996)

!.....set the electron mean molecular weight
      memw = 1.e0 / Ye

!.....set the fitting constants according to temperature
      if(T.lt.8.617e-1) then
         a0 =  6.002e+19
         a1 =  2.084e+20
         a2 =  1.872e+21
         b1 =  9.383e-1
         b2 = -4.141e-1
         b3 =  5.829e-2
         cc =  5.592e+0
      else
         a0 =  6.002e+19
         a1 =  2.084e+20
         a2 =  1.872e+21
         b1 =  1.2383e+0
         b2 = -0.8141e+0
         b3 =  0.0000e+0
         cc =  4.9924e+0
      endif 

!.....fitting formulas
      lambda= T/units%me
      xi=(rho/(memw*1.e9))**(1.e0/3.e0)/lambda
      qpair=(10.748e0*lambda**2.e0 + 0.3967e0*lambda**0.5e0 + 1.005e0)
     &      **(-1.e0) * ( 1.e0 + (rho/memw) * (7.692e+7*lambda**3.e0
     &      + 9.715e+6*lambda**0.5e0)**(-1.e0))**(-0.3e0)
      fpair=( a0 + a1*xi + a2*xi**2.e0 )* EXP( -cc*xi )
     &      /( xi**3.e0 + b1/lambda + b2/lambda**(2.e0) + 
     &      b3/lambda**(3.e0) )
      g= 1.e0 - 13.04e0*(lambda**2.e0)+1.335e+2*(lambda**4.e0)
     &     + 1.534e+3*(lambda**6.e0) + 9.186e+2*(lambda**8.e0)
      
!.....pref is the conversion factor ergs to MeV, because
!     Itoh fitting formula is expressed in cgs units [ergs/cm^3/s]
!     the factor 2 is due to the fact of considering mu and tau flavours
      qep = 0.5e0 *2.e0 *(((1.e0-cv)**2.e0+(1.e0-ca)**2.e0)+
     &       ((1.e0-cv)**2.e0-(1.e0-ca)**2.e0)*qpair)
     &       *g* EXP(-2.e0/lambda)*fpair

!.....finally, the pair production can be arbitrary multiplied by a 
!     constant factor, according to the comparison with Boltztran
      qep = qep * pair_num

      end subroutine pair

!=======================================================================
!
!     LEAK: diffusion energy rate
!
!=======================================================================

      subroutine diffusion(T,tdif,stat_weight,qdif,status)

!=======================================================================
! This subroutine calculates the diffusion energy rate, assuming the
! neutrinos to be an ideal fermion gas with zero chemical potential, and
! assuming to know the neutrino diffusion timescale
!=======================================================================

      implicit none

      integer, intent(out) :: status
      real, intent(in) :: T
      real, intent(in) :: stat_weight
      real, dimension(ne), intent(in) :: tdif

      real, intent(out) :: qdif

!.......................................................................
!     Input variables:
!     T           ----> temperature [MeV]
!     tdif        ----> vector of diffusion timescale for neutrinos
!     stat_weight ----> statistical weight for neutrinos:
!                       1 for electric neutrinos
!                       1 for electric anti-neutrinos
!                       4 for mu and tau neutrinos and anti-neutrinos
!
!     Output variables
!     qdif        ----> energy diffusion rate [erg/cm^3/s]
!     status      ----> 0 if qdif is ok
!                       1 if qdif=Infinity
!
!.......................................................................

      real :: temp
      real :: pref
      integer :: k

      status = 0
      pref = 4.*units%pi*stat_weight/((units%h*units%c)**3.)
      qdif = 0.
      do k=1,ne
        if (tdif(k).gt.0.) then
          temp = exp(-E(k)/T)/((1. + exp(-E(k)/T)))
          qdif = qdif + temp * (E(k)**3. * de(k)) / tdif(k)
        else
          status = 1
        endif
      enddo

      qdif = units%MeV * qdif * pref

      return      
      end subroutine diffusion

!=======================================================================
!
!     LEAK: diffusion timescale
!
!=======================================================================

      subroutine dif_time(opdep,mfp,tdif)

!=======================================================================
! This subroutine calculates the timescale of neutrino diffusion in a
! dense environment, according to the model developed by 
! Rosswog&Liebendorfer(2003), Appendix
!=======================================================================

      implicit none

      real, dimension(ne), intent(in) :: opdep,mfp
      real, dimension(ne), intent(out) :: tdif

      integer :: k
      real :: deltax
      real, parameter :: num_dif = 6.e0

      do k=1,ne
         deltax = opdep(k) * mfp (k)
         tdif(k) = num_dif * deltax * opdep(k) / units%c
      enddo

      end subroutine dif_time

!=======================================================================

      end module mu_tau_leakage_module

!***********************************************************************
