!=======================================================================
!
!     progenitor: read in initial fe-core from Heger et al.
!     http://homepages.spa.umn.edu/~alex/stellarevolution/
!     Copyright: Matthias Liebendoerfer, 30.09.2011, GNU public license.
!
!=======================================================================

      subroutine progenitor(state)

      use input_file_module
      use state_vector_module
 
      implicit none

      type(state_vector), intent(inout) :: state

!-----------------------------------------------------------------------
!
!     Reads in stellar profile from an input file and returns a
!     corresponding state vector. The staggering of zones between
!     the a-grid (zone edges) and the b-grid (zone centers) is set
!     later in grid_initialisation.f.
!
!-----------------------------------------------------------------------

      integer :: i,ip,iq,iy,np,nq,ny,status
      real :: f,tmp1,tmp2,tmp3
      real, dimension(nnp) :: q
      real, dimension(nnp) :: amp,arp,aup,bdp,btp,byp
      real, dimension(nnq) :: amq,arq,auq,bdq,btq,byq
      real, dimension(nnq,nny) :: jy

!.....initialise input..................................................
      if (input_path%initial) call input_read_path
      if (input_state%initial) call input_read_state
      np = input_state%np
      nq = input_state%nq
      state%nq = nq
      ny = input_state%ny
      state%ny = ny
      do i=1,ny
        jy(1:nq,i) = (/(i,i=0,nq-1)/)*ny + i
        state%jy(1:nq,i) = jy(1:nq,i)
      enddo
      state%t = 0.                                                  ![s]

!.....read in Woosley & Weaver 1995 model...............................
      if ((trim(input_path%progenitorfn).eq.'s15s7b2')
     &  .or.(trim(input_path%progenitorfn).eq.'s40s7b2')
     &  .or.(trim(input_path%progenitorfn).eq.'h40z002s0')
     &  .or.(trim(input_path%progenitorfn).eq.'u40z002')
     &  .or.(trim(input_path%progenitorfn).eq.'u50z002')) then
        call s15s7b2(amp,arp,aup,bdp,btp,byp)

!.....read Heger, Woosley & Weaver 2002 model...........................
      else
        open(unit=10,file=trim(input_path%model)
     &    //trim(input_path%progenitorfn),status='old')
        do i=1,2
          read(10,*)
        enddo
        do ip=1,np
          read(10,11) amp(ip),arp(ip),aup(ip),bdp(ip),btp(ip),byp(ip)
        enddo
11      format(t9,g23.16,t34,g23.16,t59,g23.16,t84,g23.16,t109,
     &    g23.16,t259,g23.16)
        close(10)
      endif

!.....interpolate from progenitor zoning to supernova zoning............
      do ip=1,np
        q(ip) = (nq-1)*float(ip-1)/float(np-1)
      enddo
      ip = 1
      do iq=1,nq-1
        do while (q(ip+1).le.iq-1)
          ip = ip+1
        enddo
        f = (iq-1-q(ip))/(q(ip+1)-q(ip))

!.....mass..............................................................
        amq(iq) = amp(ip) + f*(amp(ip+1)-amp(ip))

!.....radius............................................................
        tmp1 = arp(ip)**3
        tmp2 = arp(ip+1)**3
        tmp3 = tmp1 + f*(tmp2-tmp1)
        arq(iq) = tmp3**(1./3.)

!.....velocity..........................................................
        tmp1 = arp(ip)**2*aup(ip)
        tmp2 = arp(ip+1)**2*aup(ip+1)
        tmp3 = tmp1 + f*(tmp2-tmp1)
        auq(iq) = tmp3/arq(iq)**2

!.....density...........................................................
        tmp1 = log(bdp(ip))
        tmp2 = log(bdp(ip+1))
        tmp3 = tmp1 + f*(tmp2-tmp1)
        bdq(iq) = exp(tmp3)

!.....temperature.......................................................
        tmp1 = log(btp(ip))
        tmp2 = log(btp(ip+1))
        tmp3 = tmp1 + f*(tmp2-tmp1)
        btq(iq) = exp(tmp3)

!.....electron fraction.................................................
        byq(iq) = byp(ip) + f*(byp(ip+1)-byp(ip))
      enddo
      amq(nq) = amp(np)
      arq(nq) = arp(np)
      auq(nq) = aup(np)
      bdq(nq) = bdp(np)
      btq(nq) = btp(np)
      byq(nq) = byp(np)

!.....fill values into state vector.....................................
      do iq=1,nq
        state%y(jy(iq,1)) = amq(iq)                                 ![g]
        state%y(jy(iq,2)) = arq(iq)                                ![cm]
        if (iq.eq.1) then
          state%y(jy(iq,3)) = 0.0                                 ![cm/s]
        else
          state%y(jy(iq,3)) = auq(iq)                              ![cm/s]
        endif
        state%y(jy(iq,4)) = amq(iq)                                 ![g]
        state%y(jy(iq,5)) = bdq(iq)                            ![g/cm^3]
        state%y(jy(iq,6)) = btq(iq) * units%kB/units%MeV          ![MeV]
        state%y(jy(iq,7)) = byq(iq)                                  ![]
        state%y(jy(iq,8)) = 1.                                       ![]
      enddo

!.....set neutrino quantities...........................................
      do iy=9,ny
        do iq=1,nq
          state%y(jy(iq,iy)) = 0.
        enddo
      enddo

!.....rewrite enclosed mass as mass differences.........................
      state%dx = 0.
      call state_newdx(state)

!.....fill remaining quantities of state vector.........................
      call state_fill_zones(state,status)
 
      end subroutine progenitor

!=======================================================================
