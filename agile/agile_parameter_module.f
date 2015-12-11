!***********************************************************************
!
!     AGILE: parameter definition
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!***********************************************************************

      module agile_parameter_module

C.....basic parameters..................................................
      integer, parameter :: nnp=500    !number of progenitor zones
      integer, parameter :: nnq=150    !number of supernova zones
      integer, parameter :: nny=12     !number of unknowns per zone
      integer, parameter :: nng=5      !number of coupled grid points

C.....derived parameters................................................
      integer, parameter :: nnc=nnq*nny   !rank of matrix
      integer, parameter :: nnb=nny*(nng+1)   !width of band_matrix
      integer, parameter :: nns=nnc*nnb   !nonzeros in sparse_matrix

      end module agile_parameter_module

!***********************************************************************
