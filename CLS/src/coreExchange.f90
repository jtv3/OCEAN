! Copyright (C) 2024 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 10-2024
!
!
program coreExchange
  use ai_kinds, only : DP
  implicit none


  integer :: ZZ, nc, lc, nptot, ntot, nspin, lmin, lmax, l, dumi, npt, maxll
  integer :: kgl, kgh, l1, l2, l3, l4, m1, m2, m3, m4, mk, k, ip, jp, istart
  integer :: i, j, nu2, nu4, kmesh(3), nk, isite, iisite, nsite
  integer, allocatable :: nproj(:)
  real(DP) :: dumr, pi, su, yp( 0 : 1000 ), avecs(3,3), omega
  real(DP), allocatable :: x(:), w(:), gk(:,:,:)
  complex(DP) :: f1, f2, f3
  complex(DP), allocatable :: cks(:,:,:), denMat(:,:), dmat(:,:)
  logical, parameter :: yes = .true.
  logical, parameter :: no = .false.
  character(len=24) :: str
  character(len=10) :: add10
  character(len=2) :: el = 'N_'
    !
  include 'sphsetnx.h.f90'
  !
  include 'sphsetx.h.f90'
  ! Currently newgetlym is only programmed for lmax = 5
!  maxll = min( maxll, 5 )
  maxll = 5
  call newgetprefs( yp, maxll, nsphpt, wsph, xsph, ysph, zsph )
  pi = 4.0d0 * atan( 1.0d0 )
  !
  open( unit=99, file='Pquadrature', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) npt
  allocate( x( npt ), w( npt ) )
  su = 0
  do i = 1, npt
     read ( 99, * ) x( i ), w( i )
     su = su + w( i )
  end do
  close( unit=99 )
  w = w * 2 / su
  
  open(unit=99,file='kmesh.ipt',status='old')
  read(99,*) kmesh(:)
  close(99)
  nk = product(kmesh)
!  write(6,*) nk

  open(unit=99,file='avecsinbohr.ipt',form='formatted',status='old')
  read(99,*) avecs(:,:)
  close(99)
  call getomega( avecs, omega )


!  open(unit=99,file='ZNL',form='formatted',status='old')
!  read(99,*) ZZ, nc, lc
!  close(99)
!  write(add10, '(A1,I3.3,A1,I2.2,A1,I2.2)') 'z', ZZ, 'n', nc, 'l', lc

  ! Later maybe add brange, etc, for consistency check
  ! and if we have metals we'll need to figure out occ/unocc
  
!  write(str,'(A8,I3.3)') 'prjfilez', ZZ
!  open(unit=99,file=str,form='formatted',status='old')
!  read(99,*) lmin, lmax, dumi, dumr
!  allocate( nproj(lmin:lmax) )
!  do l=lmin, lmax
!    read(99,*) nproj(l)
!  enddo
!  close(99)


  write(6,'(A2,A5,X,A14,A14,A14)') '##', 'N L site', 'Real (eV)', 'Imag (eV)', 'Den Trace'
  ! 
  !
!  open(unit=99,file='edgelist',form='formatted',status='old')
!  read(99,*) ZZ, nc, lc
!  close(99)

  open(unit=98,file='exx.inp',form='formatted',status='old')
!  open(unit=98,file='sitelist',form='formatted',status='old') 
  read(98,*) nsite
  do isite = 1,nsite
!    read(98,*) el, ZZ, iisite
    read(98,*) el, ZZ, nc, lc, iisite

    write(str,'(A8,I3.3)') 'prjfilez', ZZ
    write(add10, '(A1,I3.3,A1,I2.2,A1,I2.2)') 'z', ZZ, 'n', nc, 'l', lc
    
    open(unit=99,file=str,form='formatted',status='old')
    read(99,*) lmin, lmax, dumi, dumr
    allocate( nproj(lmin:lmax) )
    do l=lmin, lmax
      read(99,*) nproj(l)
    enddo
    close(99)


  write(str,'(A8,A2,I4.4)') 'parcksv.', el, iisite
  open(unit=99,file=str,form='unformatted',access='stream',status='old')
  read(99) nptot, ntot, nspin
!  write(6,*) nptot, ntot, nspin
  allocate(cks(nptot, ntot, nspin) )
  read(99) cks
  close(99)

  allocate(dmat(nptot,nptot))
  dmat = 0.0_DP
  ip = 0
  jp = 0
!  lmax = 0
!  lmin = 0
  do l = lmin, lmax
    kgl = abs(l-lc)
    kgh = l+lc
    if( kgh .ge. kgl ) then
      allocate( gk( nproj(l), nproj(l), kgl : kgh ) )
      gk = 0.0_DP
      do k = kgl, kgh, 2
        write(str, '(1a2,3i1,1a10)' ) 'gk', lc, l, k, add10
!        write(6,*) str
        open( unit=99, file=str, form='formatted', status='old' )
        rewind 99
        read ( 99, * ) gk( :, :, k )
        close(99)
      enddo
      l1 = lc; m1 = 0
      l2 = l
      l3 = l
      l4 = lc; m4 =0
      do m2 = -l, l
        do nu2 = 1, nproj(l)
          jp = jp + 1
          do m3 = -l, l
            mk = m1 - m3
            istart = ip + (m3+l)
            if ( m1 + m2 .eq. m3 + m4 ) then
              do k = kgl, kgh, 2
                if ( abs( mk ) .le. k ) then
                  call threey( l1, m1, k, mk, l3, m3, no, npt, x, w, yp, f1 )
                  call threey( l2, m2, k, mk, l4, m4, yes, npt, x, w, yp, f2 )
!                  write(6,*) l1, m1, k, mk, l3, m3, f1
!                  write(6,*) l2, m2, k, mk, l4, m4, f2
                  do nu4 = 1, nproj(l)
                    dmat(istart + nu4,jp) = dmat(istart + nu4,jp) &
                                          + gk(nu2,nu4,k) * f1 * f2 * ( 4 * pi / ( 2 * k + 1 ) )
                  enddo
                endif
              enddo
            endif
          enddo
        enddo
      enddo
      deallocate(gk)
    else
      jp = jp + (2*l+1) * nproj(l)
    endif
    ip = ip + (2*l+1) * nproj(l)
  enddo

  f3 = 0.0_DP
  do i = 1, ntot
    do j = 1, nptot
      do k = 1, nptot
        ! Exchange has a minus sign
        f3 = f3 - dmat(k,j) * cks(k,i,1) * conjg(cks(j,i,1))
      enddo
    enddo
  enddo

  allocate(denMat(nptot, nptot) )
  denMat(:,:) = 0.0_DP
  su = 0.0_DP
  do i = 1, ntot
    do j = 1, nptot
      do k = 1, nptot
        denMat(k,j) = denMat(k,j) + cks(k,i,1) * conjg(cks(j,i,1)) / real(nk,DP)
      enddo
      su = su+ (cks(j,i,1) * conjg(cks(j,i,1))) / real(nk,DP) / omega
    enddo
  enddo
!  write(6,*) 'Local den matrix trace', su
  write(6,'(A2,1X,I1.1,1X,I1.1,1X,I4.4,1X,F14.6,E14.6,F14.6)') el, nc, lc, iisite, & 
                                                       real(f3)/real(nk,DP)/omega, &
                                                       aimag(f3)/real(nk,DP)/omega, su

!  do j = 1, 5
!    write(6,*) real(denMat(1:3,j))
!  enddo

  deallocate( dmat, cks, nproj, denMat )
  enddo


end program coreExchange
