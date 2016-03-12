      module ec_module

      use egroup_module, only: ne

      implicit none

!.....array size for ec-table...........................................
      integer, parameter :: itmax = 21
      integer, parameter :: iymax = 18
      integer, parameter :: idmax = 46
      integer, parameter :: iemax = 201

!.....ec-table..........................................................
      real, dimension(idmax,iymax,itmax,ne) :: ec
      real, dimension(itmax)                :: ec_t
      real, dimension(iymax)                :: ec_ye
      real, dimension(idmax)                :: ec_d

!.....ec-table energy grid..............................................
      real, parameter :: ec_erbmin=0.
      real, parameter :: ec_erbmax=100.

      end module
